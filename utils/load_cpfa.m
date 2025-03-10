function img_cropped = load_cpfa(filename)
    rng(0);
    img = im2double(imread(filename));
    img_size = size(img);
    h = 1024;
    w = 1024; % size for visualization
    height = img_size(1);
    width = img_size(2);

    if height <= h && width <= w
        img_cropped = img;
    else
        start_row = randi([1, height-h+1]);
        start_col = randi([1, width-w+1]);
        img_cropped = img(start_row:start_row+h-1, start_col:start_col+w-1, :);
    end
end
