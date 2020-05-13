all:
	gcc -std=c99 -lm knnd_test.c knnd.c vec.c -o knnd_test
	gcc -std=c99 -lm -lz run_approx_knn.c knnd.c vec.c -o run_approx_knn
