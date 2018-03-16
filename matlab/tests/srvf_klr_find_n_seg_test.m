%% Test 1 - Simple grid
w = [1, 0, 1;
     0, 0, 0;
     1, 0, 1];
w = flipud(w');

seg = srvf_klr_find_n_seg(2, 2, w, 0);
cor = [3, 3];

segs = sortrows(seg);
cors = sortrows(cor);
assert( all(segs(:) == cors(:)) );

%% Test 2 - Multiple segments
w = [0, 1, 0, 0, 0;
     0, 0, 1, 0, 0;
     0, 0, 0, 0, 1;
     0, 0, 0, 0, 0;
     1, 0, 0, 0, 1];
w = flipud(w)';

seg = srvf_klr_find_n_seg(2, 2, w, 0);
cor = [2, 5; 3, 4; 5, 3];

segs = sortrows(seg);
cors = sortrows(cor);
assert( all(segs(:) == cors(:)) );

%% Test 3 - Purely horizontal
w = [0 0 0 0;
     1 0 0 1;
     1 0 0 0];
w = flipud(w)';

seg = srvf_klr_find_n_seg(2, 2, w, 0);
cor = [4, 2];

segs = sortrows(seg);
cors = sortrows(cor);
assert( all(segs(:) == cors(:)) );

%% Test 4 - Empty grid
w = [0 0 0 0;
     0 0 0 0;
     0 0 0 0];
w = flipud(w)';

seg = srvf_klr_find_n_seg(1, 1, w, 0);
cor = [1, 4; 
       2, 4; 
       3, 4; 
       4, 4; 
       5, 4; 
       5, 3; 
       5, 2;
       5, 1];

segs = sortrows(seg);
cors = sortrows(cor);
assert( all(segs(:) == cors(:)) );