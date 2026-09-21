import argparse
import os
import shutil
import subprocess
import sys
import time
import threading
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor

def _normalize_rel(path_obj):
    return str(path_obj).replace("\\", "/")

def _collect_all_images_sorted(image_root):
    """
    Scan all images and sort them globally:
    1. Sort by parent-directory path first (so directories are grouped).
    2. Then sort by filename.
    """
    images_abs = [p for p in image_root.rglob("*.png") if p.is_file()]
    # Convert to relative paths.
    images_rel = [_normalize_rel(p.relative_to(image_root)) for p in images_abs]
    
    # Primary sorting logic.
    images_rel.sort(key=lambda x: (str(Path(x).parent), x))
    
    return images_rel

class TaskManager:
    """Thread-safe task manager using a producer-consumer model."""
    def __init__(self, all_images, chunk_size):
        self.all_images = all_images
        self.total = len(all_images)
        self.chunk_size = chunk_size
        self.current_idx = 0
        self.lock = threading.Lock()
        
    def get_next_batch(self):
        with self.lock:
            if self.current_idx >= self.total:
                return None # No tasks remain.
            
            # Retrieve a task batch.
            start = self.current_idx
            end = min(self.current_idx + self.chunk_size, self.total)
            batch = self.all_images[start:end]
            
            # Advance the task pointer.
            self.current_idx = end
            return batch, start, end

def copy_file_task(src, dst):
    try:
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
        return 1
    except Exception as e:
        print(f"[Merge Error] {src} -> {dst}: {e}")
        return 0

def worker_thread(gpu_id, task_manager, image_root, save_dir, workdir, project_root, use_cuda):
    """
    GPU worker thread: continuously acquire tasks from the manager until completion.
    """
    processed_count = 0
    
    # Create a separate temporary directory for each GPU.
    temp_dir = Path(save_dir) / f"gpu_{gpu_id}"
    temp_dir.mkdir(parents=True, exist_ok=True)
    
    while True:
        # 1. Acquire a task batch (atomic operation).
        task_data = task_manager.get_next_batch()
        if task_data is None:
            break # All tasks are complete.
            
        batch_imgs, start_idx, end_idx = task_data
        
        # 2. Pre-create subdirectories to avoid save failures in demo_real_save.py.
        for img_rel in batch_imgs:
            subdir = Path(img_rel).parent
            if str(subdir) != ".":
                (temp_dir / subdir).mkdir(parents=True, exist_ok=True)
        
        # 3. Build the command that invokes the subprocess.
        cmd = [
            sys.executable,
            str(workdir / "demo_real_save.py"),
            "--image-root", str(image_root),
            "--image-names", ",".join(batch_imgs),
            "--gpu-id", str(gpu_id),
            "--save-dir", str(temp_dir),
        ]
        if use_cuda:
            cmd.append("--use-cuda")
            
        env = os.environ.copy()
        env["USE_EXPLICIT_DFT"] = "1"
        env["USE_MATLAB_ENGINE_IFFT"] = "0"
        
        try:
            # Wait for this batch; then immediately acquire the next batch.
            ret = subprocess.call(cmd, cwd=str(project_root), env=env)
            if ret != 0:
                print(f"[Error] GPU {gpu_id} failed batch starting with {batch_imgs[0]}")
            else:
                processed_count += len(batch_imgs)
        except Exception as e:
            print(f"[Exception] GPU {gpu_id}: {e}")
            
    return processed_count, temp_dir

def main():
    parser = argparse.ArgumentParser(description="Global Stream Multi-GPU Launcher (Smart Balanced)")
    parser.add_argument("--gpus", type=str, default="0", help="gpu ids")
    parser.add_argument("--image-root", type=str, default=str(Path("data") / "real"))
    parser.add_argument("--save-dir", type=str, default="real_results_stream")
    parser.add_argument("--use-cuda", action="store_true")
    # The default upper limit is 10 and is reduced automatically for small inputs.
    parser.add_argument("--cmd-chunk-size", type=int, default=10, help="max images per subprocess call") 
    args = parser.parse_args()

    workdir = Path(__file__).resolve().parent
    project_root = workdir
    
    image_root = Path(args.image_root)
    if not image_root.is_absolute():
        image_root = (project_root / image_root).resolve()
        
    gpu_ids = [g.strip() for g in args.gpus.split(",") if g.strip()]
    if not gpu_ids: raise RuntimeError("No GPUs")
    
    # 1. Global scan and sort.
    print("Scanning images...")
    all_images = _collect_all_images_sorted(image_root)
    total_imgs = len(all_images)
    num_gpus = len(gpu_ids)
    print(f"Total images: {total_imgs} | Available GPUs: {num_gpus} ({gpu_ids})")
    
    if total_imgs == 0: return

    # =========================================================================
    # 2. Adaptive chunking strategy.
    # This prevents the first GPU from consuming all tasks for small inputs and
    # preserves streaming load balancing.
    # =========================================================================
    if total_imgs <= num_gpus * 2:
        # Very small input: force a chunk size of 1 so each GPU can receive work.
        dynamic_chunk = 1
    else:
        # Typical input: choose a size that lets all GPUs run several rounds.
        # The value 4 is empirical and targets at least four task acquisitions per GPU.
        ideal_chunk = max(1, total_imgs // (num_gpus * 4))
        # Never exceed the user-specified limit (10 by default).
        dynamic_chunk = min(args.cmd_chunk_size, ideal_chunk)
        
    print(f"Strategy: Chunk Size set to [{dynamic_chunk}].")
    print(f"Details:  Each GPU will take {dynamic_chunk} image(s), process them, then come back for more.")
    
    # Initialize the task manager.
    manager = TaskManager(all_images, dynamic_chunk)
    
    t0 = time.perf_counter()
    
    # 3. Start multithreaded workers (the actual parallel entry point).
    worker_results = []
    # Use a thread pool with one worker per GPU.
    with ThreadPoolExecutor(max_workers=num_gpus) as executor:
        futures = []
        for gid in gpu_ids:
            # Bind each thread to one GPU ID and run worker_thread independently.
            f = executor.submit(
                worker_thread,
                gid, manager, image_root, args.save_dir, workdir, project_root, args.use_cuda
            )
            futures.append(f)
            
        # Wait for all threads to complete.
        for f in futures:
            worker_results.append(f.result())

    # 4. Merge results.
    print("="*60)
    print("Merging results...")
    final_dir = project_root / args.save_dir
    final_dir.mkdir(parents=True, exist_ok=True)
    
    temp_dirs_to_clean = []
    copy_tasks = []
    
    for count, temp_dir in worker_results:
        temp_dirs_to_clean.append(temp_dir)
        if not temp_dir.exists(): continue
        for f in temp_dir.rglob("*"):
            if f.is_file():
                rel = f.relative_to(temp_dir)
                dst = final_dir / rel
                copy_tasks.append((f, dst))
    
    if copy_tasks:
        # Copy in parallel to improve I/O throughput.
        with ThreadPoolExecutor(max_workers=16) as pool:
            results = [pool.submit(copy_file_task, s, d) for s, d in copy_tasks]
            _ = [r.result() for r in results]
    
    # Clean up temporary directories.
    for td in temp_dirs_to_clean:
        # Attempt to remove the parent gpu_x directory.
        shutil.rmtree(td, ignore_errors=True)

    print(f"All Completed. Total Time: {time.perf_counter() - t0:.2f}s")
    print(f"Output Directory: {final_dir}")

if __name__ == "__main__":
    main()
