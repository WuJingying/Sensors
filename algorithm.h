/** file algorithm.h */


#ifndef ALGORITHM_H_
#define ALGORITHM_H_

#define true 1
#define false 0
#define FS 25
#define BUFFER_SIZE  (FS* 12) 
#define HR_FIFO_SIZE 7
#define MA4_SIZE  4 // DO NOT CHANGE
#define HAMMING_SIZE  5// DO NOT CHANGE
#define min(x,y) ((x) < (y) ? (x) : (y))

const double auw_hamm[HAMMING_SIZE] = { 40.9600000000000,276.480000000000,512,276.480000000000,40.9600000000000 }; 
//Hamm=long16(512* hamming(5)');

//uch_spo2_table is computed as  -45.060*ratioAverage* ratioAverage + 30.354 *ratioAverage + 94.845 ;
// 存放准确的血氧值
const int uch_spo2_table[184] = { 95, 95, 95, 96, 96, 96, 97, 97, 97, 97, 97, 98, 98, 98, 98, 98, 99, 99, 99, 99,
	99, 99, 99, 99, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100,
	100, 100, 100, 100, 99, 99, 99, 99, 99, 99, 99, 99, 98, 98, 98, 98, 98, 98, 97, 97,
	97, 97, 96, 96, 96, 96, 95, 95, 95, 94, 94, 94, 93, 93, 93, 92, 92, 92, 91, 91,
	90, 90, 89, 89, 89, 88, 88, 87, 87, 86, 86, 85, 85, 84, 84, 83, 82, 82, 81, 81,
	80, 80, 79, 78, 78, 77, 76, 76, 75, 74, 74, 73, 72, 72, 71, 70, 69, 69, 68, 67,
	66, 66, 65, 64, 63, 62, 62, 61, 60, 59, 58, 57, 56, 56, 55, 54, 53, 52, 51, 50,
	49, 48, 47, 46, 45, 44, 43, 42, 41, 40, 39, 38, 37, 36, 35, 34, 33, 31, 30, 29,
	28, 27, 26, 25, 23, 22, 21, 20, 19, 17, 16, 15, 14, 12, 11, 10, 9, 7, 6, 5,
	3, 2, 1 }; 


//计算心率&血氧
void maxim_heart_rate_and_oxygen_saturation(double* pun_ir_buffer, int n_ir_buffer_length, double* pun_red_buffer, int* pn_spo2, int* pch_spo2_valid, int* pn_heart_rate, int* pch_hr_valid);

//找出距离大于等于MIN_DISTANCE的、高于MIN_HEIGHT的最多MAX_NUM个波峰
void maxim_find_peaks(int* pn_locs, int* pn_npks, double* pn_x, int n_size, double n_min_height, int  n_min_distance, int n_max_num);

//寻找大于min_height的波峰
void maxim_peaks_above_min_height(int* pn_locs, int* pn_npks, double* pn_x, int n_size, double n_min_height);

//去除间隔小于MIN_DISTANCE的波峰
void maxim_remove_close_peaks(int* pn_locs, int* pn_npks, double* pn_x, int n_min_distance);

//两个排序函数
void maxim_sort_ascend(int* pn_x, int n_size);
void maxim_sort_indices_descend(double* pn_x, int* pn_indx, int n_size);

#endif /* ALGORITHM_H_ */
