import math
import numpy as np
import pywt  #安装PyWavelets库
import argparse
import matplotlib.pyplot as plt
from scipy import signal as sig

## 对ECG采集器采到的数据滤波，定位QRS波，计算心率
## 在Terminal输入命令：python ecg_basics.py -e "163.dat"   运行程序

# 画图
def draw(data):
    plt.figure(figsize=(16, 2))
    plt.plot(data)
    plt.show()


# 过滤肌电干扰、工频干扰
# 每个人的皮肤表面都有大约30mV左右的电势存在。 这一数据一般保持在25-35mV之间，
# 最大值为5mV的干扰足以对ECG信号产生干扰，主要是由于肌肉纤维的颤抖导致体表的电位发生变化，
# 这样使得体表的电极贴片测得的电势差受到影响，造成的干扰叫做肌电干扰，持续时间较短，
# 使得ECG信号波形产生细小的波纹，这种噪声的频率分布范围比较宽广，一般在零到一万赫兹之间，
# 但是相对来说分布在30-300Hz之间的更多，频率特性相当于白噪声。
# 正因为它的频率分布范围较广，传统的消噪方法取得的结果不是很理想，
# 这是由于一般的消噪方法只能很好地去除分布在特定频率段的噪声，所以肌电干扰一直是研究的难点。
# 工频干扰是在采集阶段造成的，国内一般是指50Hz的电源线干扰和各次谐波干扰，
# MIT-BIH心率失常数据库中的数据均指频率在60Hz左右的各次谐波，这是因为美国的市电是60Hz。
# 幅值大约在0-0.4mV，根据不同的情况，一般相当于R波幅值的5%-40%，
# 这是一种一般采样到的信号都会碰到的较为常见的干扰，使得ECG信号的SNR下降，甚至淹没原始信号。
# 另外，测量时由于不对称的线路产生的50Hz的工频干扰以及周围环境的电磁干扰都对ECG信号正确采集造成干扰，有待进行软件处理。
def filter_01(data):
    coeffs = pywt.wavedec(data=data, wavelet='db5', level=9)
    cA9, cD9, cD8, cD7, cD6, cD5, cD4, cD3, cD2, cD1 = coeffs
    threshold = (np.median(np.abs(cD1)) / 0.6745) * (np.sqrt(2 * np.log(len(cD1))))
    cD1.fill(0)
    cD2.fill(0)
    for i in range(1, len(coeffs) - 2):
        coeffs[i] = pywt.threshold(coeffs[i], threshold)
    rdata = pywt.waverec(coeffs=coeffs, wavelet='db5')
    return rdata


# 过滤基线漂移
# 基线漂移属于低频干扰，呼吸的节奏、四肢动作以及前端处理电路设计都有可能造成基线漂移
# ECG信号的一般采用是粘贴式或吸球式电极来采集信号，那么存在于体表与电极之间的接触类电阻和放大器产生的
# 输入阻抗，两者会形成一个分压网络，此电压将会使ECG信号发生漂移，而人体稍微的运动便会使这个接触电阻发
# 生改变，从而对原始的ECG信号产生影响。对ECG信号的影响最大，在幅值上有明显的抖动，
# 被这类噪声干扰的ECG信号变得扭曲，上下抖动，从而难以准确定位信号中特征点的位置，
# 为了得到有效的ECG信号，这类噪声是必须消除的。
def filter_02(data, fs=500):
    data = np.array(data)
    winsize = int(round(0.2 * fs))
    if winsize % 2 == 0:
        winsize += 1
    baseline_estimate = sig.medfilt(data, kernel_size=winsize)
    winsize = int(round(0.6 * fs))
    if winsize % 2 == 0:
        winsize += 1
    baseline_estimate = sig.medfilt(baseline_estimate, kernel_size=winsize)
    ecg_blr = data - baseline_estimate
    return ecg_blr.tolist()


# QRS波定位
def getQRS(data):
    count_q = []
    count_r = []
    count_s = []

    T = 256
    N = 24
    rE = T // 3  # 整数除法
    E = T // 7
    x = data[:]
    x = x.astype("float")
    x = (x - np.mean(x)) / np.std(x)
    plt.figure(num=None, figsize=(16, 2), dpi=80)

    x1 = sig.lfilter([1, 0, 0, 0, 0, 0, -2, 0, 0, 0, 0, 0, 1], [1, -2, 1], x)
    x2 = sig.lfilter(
        [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1], [1, 1],
        x1)
    plt.plot(x2[0:3000])

    x3 = np.zeros(x.shape)

    for i in range(2, len(x2) - 2):
        x3[i] = (-1 * x2[i - 2] - 2 * x2[i - 1] + 2 * x2[i + 1] + x2[i + 2]) / (8 * T)

    x4 = x3 * x3
    x5 = np.zeros(x.shape)
    for i in range(N, len(x4) - N):
        for j in range(N):
            x5[i] += x4[i - j]
    x5 = x5 / N

    peaki = x5[0]
    spki = 0
    npki = 0
    c = 0
    peak = [0]
    threshold1 = spki
    pk = []
    for i in range(1, len(x5)):
        if x5[i] > peaki:
            peaki = x5[i]

        npki = ((npki * (i - 1)) + x5[i]) / i
        spki = ((spki * (i - 1)) + x5[i]) / i
        spki = 0.875 * spki + 0.125 * peaki
        npki = 0.875 * npki + 0.125 * peaki

        threshold1 = npki + 0.25 * (spki - npki)
        threshold2 = 0.5 * threshold1

        if (x5[i] >= threshold2):

            if (peak[-1] + N < i):
                peak.append(i)
                pk.append(x5[i])

    p = np.zeros(len(x5))
    rPeak = []
    Q = np.zeros(2)
    S = np.zeros(2)
    THR = 50
    for i in peak:
        if (i == 0 or i < 2 * rE):
            continue
        p[i] = 1

        ind = np.argmax(x2[i - rE:i + rE])
        maxIndexR = (ind + i - rE)
        rPeak.append(maxIndexR)
        plt.plot(maxIndexR, x2[maxIndexR], 'ro', markersize=12)
        count_r.append(maxIndexR)
        prevDiffQ = 0
        prevDiffS = 0
        #  找Q点
        for i in range(1, THR):
            Q[0] = x2[maxIndexR - i]
            Q[1] = x2[maxIndexR - (i + 1)]

            diffQ = Q[0] - Q[1]

            if (diffQ < prevDiffQ):
                minIndexQ = maxIndexR - i
                break
            prevDiffQ = diffQ / 5
        plt.plot(minIndexQ, x2[minIndexQ], 'bo', markersize=6)
        count_q.append(minIndexQ)

        #  找S点
        for i in range(1, THR):

            S[0] = x2[maxIndexR + i]
            S[1] = x2[maxIndexR + (i + 1)]

            diffS = S[0] - S[1]

            if (diffS < prevDiffS):
                minIndexS = maxIndexR + i
                break
            prevDiffS = diffS / 5
        plt.plot(minIndexS, x2[minIndexS], 'go', markersize=6)
        count_s.append(minIndexS)
    rPeak = np.unique(rPeak)

    plt.xlabel('time')
    plt.show()
    return count_q, count_r, count_s


# 按12位取数据
def read_uint12(data_chunk):
    data = np.frombuffer(data_chunk, dtype=np.uint8)
    fst_uint8, mid_uint8, lst_uint8 = np.reshape(data, (data.shape[0] // 3, 3)).astype(np.uint16).T
    fst_uint12 = (fst_uint8 << 4) + (mid_uint8 >> 4)
    snd_uint12 = ((mid_uint8 % 16) << 8) + lst_uint8
    return np.reshape(np.concatenate((fst_uint12[:, None], snd_uint12[:, None]), axis=1), 2 * fst_uint12.shape[0])


# 取中段time_s秒数据
def getdata_10(data, sampling_rate):
    time_s = 30  # 可以改
    rdata = data
    signal_len = len(data)
    signal_len_10 = sampling_rate * time_s
    if signal_len < signal_len_10:
        pass
    else:
        offset = math.floor((signal_len - signal_len_10) / 2)
        rdata = data[offset:signal_len_10 + offset]
    return rdata


# 获取平均心率
def get_heart_rate_mean(ecg_data_path):
    bit_len = 12  # 单个信号位长度，由数据集决定
    sampling_rate = 360  # 采样率 ，可以修改

    f = open(ecg_data_path, 'rb')
    signal_bit = f.read()
    f.close()

    signal = read_uint12(signal_bit)  # 用np读数据，做转换 8位 to 12位
    signal_abc = np.array(signal).reshape(math.floor(len(signal) / 3), 3)  # 分开3个通道
    signal_10 = getdata_10(signal_abc[:, 1], sampling_rate)  # 取2号通道10秒
    f0_signal_b = filter_01(signal_10)  # 拆分导联，过滤肌电干扰、工频干扰
    f1_signal_b = filter_02(f0_signal_b, sampling_rate)  # 过滤基线漂移
    f1_signal_b = np.array(f1_signal_b)
    # draw(f1_signal_b)  # draw一下看看吧
    count_q, count_r, count_s = getQRS(f1_signal_b)  # 定位QRS
    # 计算心率
    seconds = len(f1_signal_b) / sampling_rate
    heart_rate_mean = len(count_r) * 60 / seconds
    print(heart_rate_mean)  # 输出平均心率
    return heart_rate_mean


## 主程序
if __name__ == '__main__':
    ecg_data_path = "163.dat"
    ap = argparse.ArgumentParser()
    ap.add_argument("-e", "--ecg", required=True, help="path to the input ecg data")
    args = vars(ap.parse_args())
    ecg_data_path = args['ecg']
    get_heart_rate_mean(ecg_data_path)