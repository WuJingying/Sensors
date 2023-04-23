/** file algorithm.cpp */

#include "algorithm.h"

// ����&Ѫ�����㺯��
void maxim_heart_rate_and_oxygen_saturation(double* pun_ir_buffer, int n_ir_buffer_length, double* pun_red_buffer, int* pn_spo2, int* pch_spo2_valid, int* pn_heart_rate, int* pch_hr_valid)
/*
* ͨ�����PPG���ڵķ�ֵ����Ӧ�ĺ��/������źŵĽ���/ֱ�����֣������SPO2�ı��ʡ�
* ���㷨��Ŀ�����ֱ�M0/M3�����ڼĴ��������SPO2�Ĺ�ʽû�дﵽ��ȷ�ȡ�
* ��ˣ�Ԥ�Ƚ�׼ȷ��SPO2������uch_spo2_table[]�У��Թ��Ƚϡ�
*/
{
    double un_ir_mean, an_x[BUFFER_SIZE], an_y[BUFFER_SIZE], un_only_once;

    //ȥ���ź��е�ֱ��
    un_ir_mean = 0;
    for (int k = 0; k < n_ir_buffer_length; k++) un_ir_mean += pun_ir_buffer[k]; //����洢ir���ܺ�
    un_ir_mean = un_ir_mean / n_ir_buffer_length;  //un_ir_mean��ʾir��ƽ��ֵ

    for (int k = 0; k < n_ir_buffer_length; k++)  an_x[k] = pun_ir_buffer[k] - un_ir_mean; //ÿ��Ԫ����ƽ��ֵ�Ĳ�

    // 4����ƽ������ֵ��ƽ����
    double n_denom = 0.0;
    for (int k = 0; k < BUFFER_SIZE - MA4_SIZE; k++) {
        n_denom = (an_x[k] + an_x[k + 1] + an_x[k + 2] + an_x[k + 3]);
        an_x[k] = n_denom / 4;
    } // ÿ4����һ��ƽ��ֵ�����ƽ����ĺ����źŵĲ�ֵ


    //��ֵ������Ǽ�������֮��Ĳ����ɲ�����ʹ�䲨�����߸�ƽ��
    // һ�ײ��
    double an_dx[BUFFER_SIZE];
    for (int k = 0; k < BUFFER_SIZE - MA4_SIZE - 1; k++)
        an_dx[k] = (an_x[k + 1] - an_x[k]);
    // ���ײ��
    for (int k = 0; k < BUFFER_SIZE - MA4_SIZE - 2; k++) {
        an_dx[k] = (an_dx[k] + an_dx[k + 1]) / 2;
    }

    /*�������źŴ�������У�ÿ��FFT�任ֻ�ܶ����޳��ȵ�ʱ�����ݽ��б任��
    ��ˣ���Ҫ��ʱ���źŽ����źŽضϡ���ʹ�������źţ�����ضϵ�ʱ�䳤�Ȳ������ڵ������������ڽضϣ�
    ��ô����ȡ����źŽ������й©��Ϊ�˽����й©�����ٵ���С�̶ȣ�������Ҫʹ�ü�Ȩ������Ҳ�д�������
    �Ӵ���Ҫ��Ϊ��ʹʱ���ź��ƺ����õ�����FFT�����������Ҫ�󣬼���й©��*/

    // hamming window����һ�����������������ĳһ�����з���ֵ���������������Ϊ0
    //��ת���Σ��������ò���̽����̽�⵽����
    for (int i = 0; i < BUFFER_SIZE - HAMMING_SIZE - MA4_SIZE - 2; i++) {
        double s = 0;
        for (int k = i; k < i + HAMMING_SIZE; k++) {
            s -= an_dx[k] * auw_hamm[k - i]; //an_dx[]�洢������ir����
        }
        an_dx[i] = s / 1146; //1146��auw_hamm����Ԫ��֮��
    }


    double n_th1 = 0; // ��ֵ����
    for (int k = 0; k < BUFFER_SIZE - HAMMING_SIZE; k++) {
        n_th1 += ((an_dx[k] > 0) ? an_dx[k] : (0 - an_dx[k]));
    }
    n_th1 = n_th1 / (BUFFER_SIZE - HAMMING_SIZE);
    int an_dx_peak_locs[15]; //����λ��
    int n_npks;

    //�ҳ�������ڵ���8�ġ�����n_th1�����5������
    maxim_find_peaks(an_dx_peak_locs, &n_npks, an_dx, BUFFER_SIZE - HAMMING_SIZE, n_th1, 8, 5);
    //an_dx_peak_locs������ҵ������λ�ã�
    //n_npks���ҵ������������
    //an_dx��������ir���ݣ�
    //BUFFER_SIZE-HAMMING_SIZE��һ��ȡ��Ƭ�εĴ�С��

    int n_peak_interval_sum = 0;
    if (n_npks >= 2) { //�����������2������Լ�������
        for (int k = 1; k < n_npks; k++)
            n_peak_interval_sum += (an_dx_peak_locs[k] - an_dx_peak_locs[k - 1]); //���֮����֮��
        n_peak_interval_sum = n_peak_interval_sum / (n_npks - 1);  //һ��PPG�ź��еķ��֮��ƽ���������
        // �������ʹ�ʽ
        *pn_heart_rate = int(FS * 60 / n_peak_interval_sum);// һ���Ӳɼ�6000����������������֮��ƽ��������
        *pch_hr_valid = 1; //���ʼ�����Ч��־��1
    }
    else { //��������2���޷���������
        *pn_heart_rate = -999;
        *pch_hr_valid = 0; //���ʼ�����Ч��־��0
    }
    //���ϣ����ʼ������


    // ��ʼ���ݲ���λ������
    int an_ir_valley_locs[15]; //���ir���ݵĲ���λ�õĳ�ʼ��
    for (int k = 0; k < n_npks; k++) {
        an_ir_valley_locs[k] = an_dx_peak_locs[k] + int(HAMMING_SIZE / 2);
    }


    // ԭʼ����: RED(=y) and IR(=X)
    // ͨ��IR��RED�����DC��AC
    for (int k = 0; k < n_ir_buffer_length; k++) {
        an_x[k] = pun_ir_buffer[k]; //ir����֮��
        an_y[k] = pun_red_buffer[k]; //red����֮��
    }

    // ��ȷ����λ�� �Լ�Сspo2�������
    int n_exact_ir_valley_locs_count = 0; //ir����׼ȷ�Ĳ�������
    int m = 0;
    double n_c_min = 0.0;
    int an_exact_ir_valley_locs[15]; //ir����׼ȷ�Ĳ���λ��
    for (int k = 0; k < n_npks; k++) {
        un_only_once = 1;
        m = an_ir_valley_locs[k]; //mΪir���ݲ���λ�ó�ʼֵ
        n_c_min = 16777216;//2^24; //ĳ����ֵ
        if (m + 5 < BUFFER_SIZE - HAMMING_SIZE && m - 5 > 0) {
            for (int i = m - 5; i < m + 5; i++) { //����mλ��Ϊ���ģ�5Ϊ�뾶�Ĵ����У�Ѱ��irֵ��С��λ��
                if (an_x[i] < n_c_min) {
                    if (un_only_once > 0) {
                        un_only_once = 0;
                    }
                    n_c_min = an_x[i];
                    an_exact_ir_valley_locs[k] = i; //��¼׼ȷ�Ĳ���λ��
                }
            }
            if (un_only_once == 0) n_exact_ir_valley_locs_count++; //׼ȷ�ҵ���һ������
        }
    }
    // ����С��2�������޷�����spo2���������˽���
    if (n_exact_ir_valley_locs_count < 2) {
        *pn_spo2 = -999;
        *pch_spo2_valid = 0; //Ѫ��������Ч��־��0
        return;
    }

    // ���ȡ�2������������
    // ��֣�ʹred��ir���ݱ��ƽ��
    for (int k = 0; k < BUFFER_SIZE - MA4_SIZE; k++) {
        an_x[k] = (an_x[k] + an_x[k + 1] + an_x[k + 2] + an_x[k + 3]) / 4; //ir����
        an_y[k] = (an_y[k] + an_y[k + 1] + an_y[k + 2] + an_y[k + 3]) / 4; //red����
    }


    //ʹ��ir����׼ȷ����λ�ã�����ֱ��DC�ͽ���AC
    //Ѱ����������֮��AC/DC���ݵ����ֵ
    int n_ratio_average = 0;
    int n_i_ratio_count = 0;
    int an_ratio[5];
    for (int k = 0; k < 5; k++) an_ratio[k] = 0; //��ʼ��an_ratio����
    for (int k = 0; k < n_exact_ir_valley_locs_count; k++) {
        if (an_exact_ir_valley_locs[k] > BUFFER_SIZE) { //���ir����׼ȷ����λ��>300,��������ir���鷶Χ
            *pn_spo2 = -999; //������SPO2��Ч����������
            *pch_spo2_valid = 0;
            return;
        }
    }

    // �ֱ�Ѱ��ir���ݺ�red������������֮������ֵ 
    // xϵ�У�ir����    yϵ�У�red����
    int n_x_dc_max = 0, n_y_dc_max = 0; //��¼���ֵ
    int n_x_dc_max_idx, n_y_dc_max_idx = 0; //��¼���ֵ���±�
    int n_x_ac = 0, n_y_ac = 0;
    int n_nume = 0;
    for (int k = 0; k < n_exact_ir_valley_locs_count - 1; k++) {
        n_y_dc_max = -16777216;
        n_x_dc_max = -16777216; //��ʼ����������
        if (an_exact_ir_valley_locs[k + 1] - an_exact_ir_valley_locs[k] > 10) { //���ir���� �������ȵľ��룾10
            for (int i = an_exact_ir_valley_locs[k]; i < an_exact_ir_valley_locs[k + 1]; i++) {
                if (an_x[i] > n_x_dc_max) { n_x_dc_max = an_x[i]; n_x_dc_max_idx = i; }
                if (an_y[i] > n_y_dc_max) { n_y_dc_max = an_y[i]; n_y_dc_max_idx = i; }
            }//�Ҳ��ȵ����ֵ��ir������red���ݲ���λ��һ�£�����ֻ����ir���ݵĲ���λ��

            //����f(t)=dc(t)+ac(t) �ֱ���dc��ac
            //red����
            n_y_ac = (an_y[an_exact_ir_valley_locs[k + 1]] - an_y[an_exact_ir_valley_locs[k]]) * (n_y_dc_max_idx - an_exact_ir_valley_locs[k]);
            n_y_ac = an_y[an_exact_ir_valley_locs[k]] + n_y_ac / (an_exact_ir_valley_locs[k + 1] - an_exact_ir_valley_locs[k]);
            n_y_ac = an_y[n_y_dc_max_idx] - n_y_ac; // ��ԭ������ȥ��ֱ��dc
            //ir����
            n_x_ac = (an_x[an_exact_ir_valley_locs[k + 1]] - an_x[an_exact_ir_valley_locs[k]]) * (n_x_dc_max_idx - an_exact_ir_valley_locs[k]);
            n_x_ac = an_x[an_exact_ir_valley_locs[k]] + n_x_ac / (an_exact_ir_valley_locs[k + 1] - an_exact_ir_valley_locs[k]);
            n_x_ac = an_x[n_y_dc_max_idx] - n_x_ac; // ��ԭ������ȥ��ֱ��dc

            n_nume = (n_y_ac * n_x_dc_max) >> 7; //����Ϊfloat����
            n_denom = (n_x_ac * n_y_dc_max) >> 7;
            if (n_denom > 0 && n_i_ratio_count < 5 && n_nume != 0)
            {
                // ��ʽ��an_ratio = ( n_y_ac *n_x_dc_max) / ( n_x_ac *n_y_dc_max) ;
                an_ratio[n_i_ratio_count] = (n_nume * 100) / n_denom;
                n_i_ratio_count++;
            }
        }
    }
    //���ϣ������Ѫ����ʽSpo2 = A*R*R + B*R + c (n�׶���ʽ)�е�Rֵ����5���������an_ratio[]�У�

    double float_SPO2 = 0.0; //SPO2��ֵ

    //�ҳ�an_ratio�������λ����Ϊ���յ�Rֵ
    maxim_sort_ascend(an_ratio, n_i_ratio_count);
    int n_middle_idx = 0;
    n_middle_idx = n_i_ratio_count / 2;
    if (n_middle_idx > 1) {
        n_ratio_average = (an_ratio[n_middle_idx - 1] + an_ratio[n_middle_idx]) / 2; //ʹ����λ��
    }
    else {
        n_ratio_average = an_ratio[n_middle_idx];
    }

    int n_spo2_calc = 0;
    if (n_ratio_average > 2 && n_ratio_average < 184) { //���2<R<184��184Ϊ��ȷѪ��ֵ����Ĵ�С
        n_spo2_calc = uch_spo2_table[n_ratio_average]; //��׼ȷֵ�Ƚ�
        *pn_spo2 = n_spo2_calc; //*pn_spo2ΪѪ��ֵ������
        *pch_spo2_valid = 1; //��Ч��־Ϊ1

        //R=n_ratio_average�����빫ʽ����SPO2ֵ
        float_SPO2 = -45.060 * n_ratio_average * n_ratio_average / 10000 + 30.354 * n_ratio_average / 100 + 94.845;
        *pn_spo2 = int(float_SPO2);
    }
    else {  // ����ó���SPO2��Ч
        *pn_spo2 = -999; 
        *pch_spo2_valid = 0; //��Ч��־Ϊ0
    }
}


void maxim_find_peaks(int* pn_locs, int* pn_npks, double* pn_x, int n_size, double n_min_height, int  n_min_distance, int n_max_num)
// �ҳ�������ڵ���MIN_DISTANCE�ġ�����MIN_HEIGHT�����MAX_NUM������
{
    // ��ȥ��С��MIN_HEIGHT�Ĳ���
    maxim_peaks_above_min_height(pn_locs, pn_npks, pn_x, n_size, n_min_height);

    //��ȥ����������Ĳ���
    maxim_remove_close_peaks(pn_locs, pn_npks, pn_x, n_min_distance);

    //
    *pn_npks = min(*pn_npks, n_max_num);
}

void maxim_peaks_above_min_height(int* pn_locs, int* pn_npks, double* pn_x, int n_size, double n_min_height)
// �ҳ����и���MIN_HEIGHT�Ĳ���
{
    int i = 1, n_width = 0;
    *pn_npks = 0;

    while (i < n_size - 1) {
        if (pn_x[i] > n_min_height && pn_x[i] > pn_x[i - 1]) { //Ѱ��Ǳ�ڷ�����Ե
            n_width = 1;
            while (i + n_width < n_size && pn_x[i] == pn_x[i + n_width])//Ѱ��ƽ����
                n_width++;
            if (pn_x[i] > pn_x[i + n_width] && (*pn_npks) < 15) { //Ѱ�ҷ���ұ�Ե
                pn_locs[(*pn_npks)++] = i; // ƽ����ȡ���Եֵ
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
// ȥ�����С��MIN_DISTANCE�ķ�ֵ
{

    int i, j, n_old_npks, n_dist;

    // �Ӵ�С���з�ֵ
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

    // ����ֵ��������
    maxim_sort_ascend(pn_locs, *pn_npks);
}

void maxim_sort_ascend(int* pn_x, int n_size)
// ������������
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
// ��������indices
{
    int i = 0, j = 0, n_temp = 0;
    for (i = 1; i < n_size; i++) {
        n_temp = pn_indx[i];
        for (j = i; j > 0 && pn_x[n_temp] > pn_x[pn_indx[j - 1]]; j--)
            pn_indx[j] = pn_indx[j - 1];
        pn_indx[j] = n_temp;
    }
}

