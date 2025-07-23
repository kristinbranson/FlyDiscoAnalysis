function posts = RemoveSmallBouts(posts,blen)
% replaced mayank's while loop with KB's suggestion - much faster same
% numerical results. 
% first replace 1 frame 0s with 1s
posts = imclose(posts,ones(1,blen));
% second replace 1 frame 1s with 0s
posts = imopen(posts,ones(1,blen));