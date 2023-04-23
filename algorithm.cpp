/** file algorithm.cpp */

#include "algorithm.h"

// 心率&血氧计算函数
void maxim_heart_rate_and_oxygen_saturation(double* pun_ir_buffer, int n_ir_buffer_length, double* pun_red_buffer, int* pn_spo2, int* pch_spo2_valid, int* pn_heart_rate, int* pch_hr_valid)
/*
* 通过检测PPG周期的峰值和相应的红光/红外光信号的交流/直流部分，计算出SPO2的比率。
* 该算法的目标是手臂M0/M3，由于寄存器溢出，SPO2的公式没有达到精确度。
* 因此，预先将准确的SPO2保存于uch_spo2_table[]中，以供比较。
*/
{
    double un_ir_mean, an_x[BUFFER_SIZE], an_y[BUFFER_SIZE], un_only_once;

    //去除信号中的直流
    un_ir_mean = 0;
    for (int k = 0; k < n_ir_buffer_length; k++) un_ir_mean += pun_ir_buffer[k]; //计算存储ir的总和
    un_ir_mean = un_ir_mean / n_ir_buffer_length;  //un_ir_mean表示ir的平均值

    for (int k = 0; k < n_ir_buffer_length; k++)  an_x[k] = pun_ir_buffer[k] - un_ir_mean; //每个元素与平均值的差

    // 4项与平均数差值的平均数
    double n_denom = 0.0;
    for (int k = 0; k < BUFFER_SIZE - MA4_SIZE; k++) {
        n_denom = (an_x[k] + an_x[k + 1] + an_x[k + 2] + an_x[k + 3]);
        an_x[k] = n_denom / 4;
    } // 每4项求一次平均值，获得平滑后的红外信号的差值


    //差分的作用是减轻数据之间的不规律波动，使其波动曲线更平稳
    // 一阶差分
    double an_dx[BUFFER_SIZE];
    for (int k = 0; k < BUFFER_SIZE - MA4_SIZE - 1; k++)
        an_dx[k] = (an_x[k + 1] - an_x[k]);
    // 二阶差分
    for (int k = 0; k < BUFFER_SIZE - MA4_SIZE - 2; k++) {
        an_dx[k] = (an_dx[k] + an_dx[k + 1]) / 2;
    }

    /*在数字信号处理过程中，每次FFT变换只能对有限长度的时域数据进行变换，
    因此，需要对时域信号进行信号截断。即使是周期信号，如果截断的时间长度不是周期的整数倍（周期截断，
    那么，截取后的信号将会存在泄漏。为了将这个泄漏误差减少到最小程度，我们需要使用加权函数，也叫窗函数。
    加窗主要是为了使时域信号似乎更好地满足FFT处理的周期性要求，减少泄漏。*/

    // hamming window构造一个函数。这个函数在某一区间有非零值，而在其余区间皆为0
    //翻转波形，就能利用波峰探测器探测到波谷
    for (int i = 0; i < BUFFER_SIZE - HAMMING_SIZE - MA4_SIZE - 2; i++) {
        double s = 0;
        for (int k = i; k < i + HAMMING_SIZE; k++) {
            s -= an_dx[k] * auw_hamm[k - i]; //an_dx[]存储处理后的ir数据
        }
        an_dx[i] = s / 1146; //1146是auw_hamm数组元素之和
    }


    double n_th1 = 0; // 阈值计算
    for (int k = 0; k < BUFFER_SIZE - HAMMING_SIZE; k++) {
        n_th1 += ((an_dx[k] > 0) ? an_dx[k] : (0 - an_dx[k]));
    }
    n_th1 = n_th1 / (BUFFER_SIZE - HAMMING_SIZE);
    int an_dx_peak_locs[15]; //波峰位置
    int n_npks;

    //找出距离大于等于8的、高于n_th1的最多5个波峰
    maxim_find_peaks(an_dx_peak_locs, &n_npks, an_dx, BUFFER_SIZE - HAMMING_SIZE, n_th1, 8, 5);
    //an_dx_peak_locs：存放找到波峰的位置；
    //n_npks：找到波峰的数量；
    //an_dx：处理后的ir数据；
    //BUFFER_SIZE-HAMMING_SIZE：一个取样片段的大小；

    int n_peak_interval_sum = 0;
    if (n_npks >= 2) { //如果波峰数≥2，则可以计算心率
        for (int k = 1; k < n_npks; k++)
            n_peak_interval_sum += (an_dx_peak_locs[k] - an_dx_peak_locs[k - 1]); //峰峰之间间隔之和
        n_peak_interval_sum = n_peak_interval_sum / (n_npks - 1);  //一段PPG信号中的峰峰之间平均间隔点数
        // 计算心率公式
        *pn_heart_rate = int(FS * 60 / n_peak_interval_sum);// 一分钟采集6000个样本，除以两峰之间平均样本数
        *pch_hr_valid = 1; //心率计算有效标志置1
    }
    else { //波峰数＜2，无法计算心率
        *pn_heart_rate = -999;
        *pch_hr_valid = 0; //心率计算有效标志置0
    }
    //以上：心率计算完成


    // 初始数据波谷位置修正
    int an_ir_valley_locs[15]; //存放ir数据的波谷位置的初始置
    for (int k = 0; k < n_npks; k++) {
        an_ir_valley_locs[k] = an_dx_peak_locs[k] + int(HAMMING_SIZE / 2);
    }


    // 原始数据: RED(=y) and IR(=X)
    // 通过IR和RED计算出DC和AC
    for (int k = 0; k < n_ir_buffer_length; k++) {
        an_x[k] = pun_ir_buffer[k]; //ir数据之和
        an_y[k] = pun_red_buffer[k]; //red数据之和
    }

    // 精确查找位置 以减小spo2计算误差
    int n_exact_ir_valley_locs_count = 0; //ir数据准确的波谷数量
    int m = 0;
    double n_c_min = 0.0;
    int an_exact_ir_valley_locs[15]; //ir数据准确的波谷位置
    for (int k = 0; k < n_npks; k++) {
        un_only_once = 1;
        m = an_ir_valley_locs[k]; //m为ir数据波谷位置初始值
        n_c_min = 16777216;//2^24; //某个阈值
        if (m + 5 < BUFFER_SIZE - HAMMING_SIZE && m - 5 > 0) {
            for (int i = m - 5; i < m + 5; i++) { //在以m位置为中心，5为半径的窗口中，寻找ir值最小的位置
                if (an_x[i] < n_c_min) {
                    if (un_only_once > 0) {
                        un_only_once = 0;
                    }
                    n_c_min = an_x[i];
                    an_exact_ir_valley_locs[k] = i; //记录准确的波谷位置
                }
            }
            if (un_only_once == 0) n_exact_ir_valley_locs_count++; //准确找到了一个波谷
        }
    }
    // 波谷小于2个，则无法计算spo2，函数到此结束
    if (n_exact_ir_valley_locs_count < 2) {
        *pn_spo2 = -999;
        *pch_spo2_valid = 0; //血氧计算有效标志置0
        return;
    }

    // 波谷≥2个，继续计算
    // 差分，使red和ir数据变得平滑
    for (int k = 0; k < BUFFER_SIZE - MA4_SIZE; k++) {
        an_x[k] = (an_x[k] + an_x[k + 1] + an_x[k + 2] + an_x[k + 3]) / 4; //ir数据
        an_y[k] = (an_y[k] + an_y[k + 1] + an_y[k + 2] + an_y[k + 3]) / 4; //red数据
    }


    //使用ir数据准确波谷位置，计算直流DC和交流AC
    //寻找两个波谷之间AC/DC数据的最大值
    int n_ratio_average = 0;
    int n_i_ratio_count = 0;
    int an_ratio[5];
    for (int k = 0; k < 5; k++) an_ratio[k] = 0; //初始化an_ratio数组
    for (int k = 0; k < n_exact_ir_valley_locs_count; k++) {
        if (an_exact_ir_valley_locs[k] > BUFFER_SIZE) { //如果ir数据准确波谷位置>300,即超出了ir数组范围
            *pn_spo2 = -999; //则计算的SPO2无效，函数结束
            *pch_spo2_valid = 0;
            return;
        }
    }

    // 分别寻找ir数据和red数据两个波谷之间的最大值 
    // x系列：ir数据    y系列：red数据
    int n_x_dc_max = 0, n_y_dc_max = 0; //记录最大值
    int n_x_dc_max_idx, n_y_dc_max_idx = 0; //记录最大值的下标
    int n_x_ac = 0, n_y_ac = 0;
    int n_nume = 0;
    for (int k = 0; k < n_exact_ir_valley_locs_count - 1; k++) {
        n_y_dc_max = -16777216;
        n_x_dc_max = -16777216; //初始化两个变量
        if (an_exact_ir_valley_locs[k + 1] - an_exact_ir_valley_locs[k] > 10) { //如果ir数据 两个波谷的距离＞10
            for (int i = an_exact_ir_valley_locs[k]; i < an_exact_ir_valley_locs[k + 1]; i++) {
                if (an_x[i] > n_x_dc_max) { n_x_dc_max = an_x[i]; n_x_dc_max_idx = i; }
                if (an_y[i] > n_y_dc_max) { n_y_dc_max = an_y[i]; n_y_dc_max_idx = i; }
            }//找波谷的最大值，ir数据与red数据波谷位置一致，所以只需存放ir数据的波谷位置

            //利用f(t)=dc(t)+ac(t) 分别求dc和ac
            //red数据
            n_y_ac = (an_y[an_exact_ir_valley_locs[k + 1]] - an_y[an_exact_ir_valley_locs[k]]) * (n_y_dc_max_idx - an_exact_ir_valley_locs[k]);
            n_y_ac = an_y[an_exact_ir_valley_locs[k]] + n_y_ac / (an_exact_ir_valley_locs[k + 1] - an_exact_ir_valley_locs[k]);
            n_y_ac = an_y[n_y_dc_max_idx] - n_y_ac; // 从原数据中去除直流dc
            //ir数据
            n_x_ac = (an_x[an_exact_ir_valley_locs[k + 1]] - an_x[an_exact_ir_valley_locs[k]]) * (n_x_dc_max_idx - an_exact_ir_valley_locs[k]);
            n_x_ac = an_x[an_exact_ir_valley_locs[k]] + n_x_ac / (an_exact_ir_valley_locs[k + 1] - an_exact_ir_valley_locs[k]);
            n_x_ac = an_x[n_y_dc_max_idx] - n_x_ac; // 从原数据中去除直流dc

            n_nume = (n_y_ac * n_x_dc_max) >> 7; //保留为float精度
            n_denom = (n_x_ac * n_y_dc_max) >> 7;
            if (n_denom > 0 && n_i_ratio_count < 5 && n_nume != 0)
            {
                // 公式：an_ratio = ( n_y_ac *n_x_dc_max) / ( n_x_ac *n_y_dc_max) ;
                an_ratio[n_i_ratio_count] = (n_nume * 100) / n_denom;
                n_i_ratio_count++;
            }
        }
    }
    //以上，计算出血氧公式Spo2 = A*R*R + B*R + c (n阶多项式)中的R值（有5个，存放于an_ratio[]中）

    double float_SPO2 = 0.0; //SPO2的值

    //找出an_ratio数组的中位数作为最终的R值
    maxim_sort_ascend(an_ratio, n_i_ratio_count);
    int n_middle_idx = 0;
    n_middle_idx = n_i_ratio_count / 2;
    if (n_middle_idx > 1) {
        n_ratio_average = (an_ratio[n_middle_idx - 1] + an_ratio[n_middle_idx]) / 2; //使用中位数
    }
    else {
        n_ratio_average = an_ratio[n_middle_idx];
    }

    int n_spo2_calc = 0;
    if (n_ratio_average > 2 && n_ratio_average < 184) { //如果2<R<184，184为精确血氧值数组的大小
        n_spo2_calc = uch_spo2_table[n_ratio_average]; //和准确值比较
        *pn_spo2 = n_spo2_calc; //*pn_spo2为血氧值计算结果
        *pch_spo2_valid = 1; //有效标志为1

        //R=n_ratio_average，代入公式计算SPO2值
        float_SPO2 = -45.060 * n_ratio_average * n_ratio_average / 10000 + 30.354 * n_ratio_average / 100 + 94.845;
        *pn_spo2 = int(float_SPO2);
    }
    else {  // 计算得出的SPO2无效
        *pn_spo2 = -999; 
        *pch_spo2_valid = 0; //有效标志为0
    }
}


void maxim_find_peaks(int* pn_locs, int* pn_npks, double* pn_x, int n_size, double n_min_height, int  n_min_distance, int n_max_num)
// 找出间隔大于等于MIN_DISTANCE的、高于MIN_HEIGHT的最多MAX_NUM个波峰
{
    // 先去除小于MIN_HEIGHT的波峰
    maxim_peaks_above_min_height(pn_locs, pn_npks, pn_x, n_size, n_min_height);

    //再去除间隔不够的波峰
    maxim_remove_close_peaks(pn_locs, pn_npks, pn_x, n_min_distance);

    //
    *pn_npks = min(*pn_npks, n_max_num);
}

void maxim_peaks_above_min_height(int* pn_locs, int* pn_npks, double* pn_x, int n_size, double n_min_height)
// 找出所有高于MIN_HEIGHT的波峰
{
    int i = 1, n_width = 0;
    *pn_npks = 0;

    while (i < n_size - 1) {
        if (pn_x[i] > n_min_height && pn_x[i] > pn_x[i - 1]) { //寻找潜在峰的左边缘
            n_width = 1;
            while (i + n_width < n_size && pn_x[i] == pn_x[i + n_width])//寻找平顶峰
                n_width++;
            if (pn_x[i] > pn_x[i + n_width] && (*pn_npks) < 15) { //寻找峰的右边缘
                pn_locs[(*pn_npks)++] = i; // 平顶峰取左边缘值
                i += n_width + 1;
            }
            else
                i += n_width;
        }
        else
            i++;
    }
}


void maxim_remove_close_peaks(int* pn_locs, int* pn_npks, double* pn_x, int n_min_distance)
// 去除间隔小于MIN_DISTANCE的峰值
{

    int i, j, n_old_npks, n_dist;

    // 从大到小排列峰值
    maxim_sort_indices_descend(pn_x, pn_locs, *pn_npks);

    for (i = -1; i < *pn_npks; i++) {
        n_old_npks = *pn_npks;
        *pn_npks = i + 1;
        for (j = i + 1; j < n_old_npks; j++) {
            n_dist = pn_locs[j] - (i == -1 ? -1 : pn_locs[i]); // lag-zero peak of autocorr is at index -1
            if (n_dist > n_min_distance || n_dist < -n_min_distance)
                pn_locs[(*pn_npks)++] = pn_locs[j];
        }
    }

    // 按峰值降序排列
    maxim_sort_ascend(pn_locs, *pn_npks);
}

void maxim_sort_ascend(int* pn_x, int n_size)
// 升序排列数组
{
    int i, j, n_temp;
    for (i = 1; i < n_size; i++) {
        n_temp = pn_x[i];
        for (j = i; j > 0 && n_temp < pn_x[j - 1]; j--)
            pn_x[j] = pn_x[j - 1];
        pn_x[j] = n_temp;
    }
}

void maxim_sort_indices_descend(double* pn_x, int* pn_indx, int n_size)
// 降序排列indices
{
    int i = 0, j = 0, n_temp = 0;
    for (i = 1; i < n_size; i++) {
        n_temp = pn_indx[i];
        for (j = i; j > 0 && pn_x[n_temp] > pn_x[pn_indx[j - 1]]; j--)
            pn_indx[j] = pn_indx[j - 1];
        pn_indx[j] = n_temp;
    }
}

