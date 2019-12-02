# -*- coding:utf-8 -*-
import csv

import math

import pandas as pd

from pandas import Series, DataFrame

import numpy as np

from numpy import *

from datetime import datetime

from timeit import timeit

import time

import matplotlib

import matplotlib.pyplot as plt

from matplotlib.pyplot import plot, savefig

from matplotlib.finance import candlestick_ochl

from matplotlib.pylab import date2num

from functools import reduce

"""
Import Data
"""

data_origin = pd.read_csv('399001.csv')
index_origin = data_origin.index.tolist()
index_reverse = index_origin[::-1]

data_origin_ = data_origin.ix[index_reverse]
data_origin_.index = index_origin


open = data_origin['open']
high = data_origin['high']
low = data_origin['low']
close = data_origin['close']
volume = data_origin['volume']

open_ = [i for i in open]
high_ = [i for i in high]
low_ = [i for i in low]
close_ = [i for i in close]
volume_ = [i for i in volume]

"""
Data Processing
"""

# date processing: 1991-04-03 ~ 2017-03-31
date = data_origin['date']
date_s2t = [datetime.strptime(d, '%Y/%m/%d') for d in date]
date_candle = [date2num(d) for d in date_s2t]


def seg(x, l):
    """ data partition """
    x = list(x)
    # print(list[1:5])
    segment = [x[i:i + l] for i in range(len(x) - l + 1)]
    # print(segment)
    return segment


# print(data_origin.apply(lambda x: seg(x, 5), axis=0))


open_full = seg(open_, 5)
high_full = seg(high_, 5)
low_full = seg(low_, 5)
close_full = seg(close_, 5)
volume_full = seg(volume_, 5)
date_full = seg(date_candle, 5)
k_line_full_df = DataFrame({'date': date_full, 'open': open_full, 'high': high_full,
                           'low': low_full, 'close': close_full, 'volume': volume_full})
df = data_origin.apply(lambda x:seg(x,5), axis=0)
# print(data_origin.ID)

# print(df)
# print('data', data_origin.ID[:5])
# print('data', data_origin)

"""
Dynamic Time Wrapping distance
"""


def merge(a, b):
    return a + b


def distance_square(history, current):
    """ calculate the distance matrix(square) between history data and current data """
    n = len(current)
    matrix_d_value = mat([[current[i] - history[j] for j in range(n)]
                          for i in range(n)]).reshape(n, n)
    matrix_d_value = matrix_d_value
    distance = multiply(matrix_d_value, matrix_d_value)
    return distance


# print(distance_square(close_full[0], close_full[0])[1])


def distance_dtw(history_df, current_df):
    """ DTW distance for k line information(open price, high price, low price,
        close price and volume) """
    open_history = reduce(merge, history_df['open'])
    open_current = reduce(merge, current_df['open'])
    high_history = reduce(merge, history_df['high'])
    high_current = reduce(merge, current_df['high'])
    low_history = reduce(merge, history_df['low'])
    low_current = reduce(merge, current_df['low'])
    close_history = reduce(merge, history_df['close'])
    close_current = reduce(merge, current_df['close'])
    volume_history = reduce(merge, history_df['volume'])
    volume_current = reduce(merge, current_df['volume'])
    # calculate the square distance matrix
    d_open = distance_square(open_history, open_current)
    d_high = distance_square(high_history, high_current)
    d_low = distance_square(low_history, low_current)
    d_close = distance_square(close_history, close_current)
    d_volume = distance_square(volume_history, volume_current)
    distance_matrix = np.sqrt(d_open + d_high + d_low + d_close + d_volume)
    # choose the minimum path
    p_i = 0
    i = 0
    j = 0
    while i < len(distance_matrix) - 1:
        c = [[distance_matrix[i + 1, j], distance_matrix[i, j + 1], distance_matrix[i + 1, j + 1]]
             if j + 1 < len(distance_matrix) else [distance_matrix[i + 1, j]]][0]
        c_min = min(c)
        p_i += c_min
        i = [i if c.index(c_min) == 1 else i + 1][0]
        j = [j if c.index(c_min) == 0 else j + 1][0]
        # print(p_i)
    return p_i


# print(distance_dtw(k_line_full_df[:1], k_line_full_df[1:2]))


""" Euclidean distance"""


def ED(history_df, current_df):
    """ Euclidean distance"""
    # get ohlc and volume
    open_history = reduce(merge, history_df['open'])
    open_current = reduce(merge, current_df['open'])
    high_history = reduce(merge, history_df['high'])
    high_current = reduce(merge, current_df['high'])
    low_history = reduce(merge, history_df['low'])
    low_current = reduce(merge, current_df['low'])
    close_history = reduce(merge, history_df['close'])
    close_current = reduce(merge, current_df['close'])
    volume_history = reduce(merge, history_df['volume'])
    volume_current = reduce(merge, current_df['volume'])
    # matrix_history = mat(open_history + high_history + low_history + close_history + volume_history)
    # matrix_current = mat(open_current + high_current + low_current + close_current + volume_current)
    # #
    # d_matrix = (matrix_current - matrix_history)
    # d_square = multiply(d_matrix, d_matrix)
    # distance = np.sqrt(sum(d_square))
    history_df_new = DataFrame({'open': open_history, 'high': high_history, 'low': low_history,
                                'close': close_history, 'volume': volume_history})
    current_df_new = DataFrame({'open': open_current, 'high': high_current, 'low': low_current,
                                'close': close_current, 'volume': volume_current})
    square_df = (current_df_new - history_df_new)**2
    # print(square_df)
    distance = np.sqrt(sum(square_df.open) + sum(square_df.high + sum(square_df.low)) +
                       sum(square_df.close) + sum(square_df.volume))
    return distance


# print(ED(k_line_full_df[:1], k_line_full_df[:1]))


""" Similarity Preference distance  """


def horizontal_shift(history, current, phi):
    """ horizontal shift similarity"""
    # calculate the mean and variance of two time series respectively
    mean_history = mean(history)
    mean_current = mean(current)
    var_history = np.std(history, ddof=1) ** 2
    var_current = np.std(current, ddof=1) ** 2
    # similarity
    alpha = math.exp(-(mean_current - mean_history) ** 2 / (phi * (var_current + var_history)))
    return alpha


# print(horizontal_shift(close_full[0], close_full[0], 0.5))


def scope(history, current, psi):
    """ scope similarity"""
    # calculate the variance of two time series respectively
    var_history = np.std(history, ddof=1) ** 2
    var_current = np.std(current, ddof=1) ** 2
    var_max = max(var_history, var_current)
    var_min = min(var_history, var_current)
    # similarity
    beta = math.exp(-(var_max - var_min) / psi * var_min)
    return beta


# print(scope(close_full[0], close_full[0], 0.5))


def correlation(history, current):
    """ correlation coefficient"""
    s_history = Series(history)
    s_current = Series(current)
    cor = s_current.corr(s_history)
    pho = [cor if cor > 0 else 0][0]
    return pho


# print(correlation(close_full[0], close_full[0]))


def similarity_preference(history_df, current_df, w1=0.1, w2=0.1, w3=0.8, phi=0.5, psi=0.5):
    """ similarity preference:  the similarity is 1 when X is the same as Y  """
    # get ohlc and volume
    open_history = reduce(merge, history_df['open'])
    open_current = reduce(merge, current_df['open'])
    high_history = reduce(merge, history_df['high'])
    high_current = reduce(merge, current_df['high'])
    low_history = reduce(merge, history_df['low'])
    low_current = reduce(merge, current_df['low'])
    close_history = reduce(merge, history_df['close'])
    close_current = reduce(merge, current_df['close'])
    volume_history = reduce(merge, history_df['volume'])
    volume_current = reduce(merge, current_df['volume'])
    #
    f_open = w1 * horizontal_shift(open_history, open_current, phi) + \
             w2 * scope(open_history, open_current, psi) + \
             w3 * correlation(open_history, open_current)
    f_high = w1 * horizontal_shift(high_history, high_current, phi) + \
             w2 * scope(high_history, high_current, psi) + \
             w3 * correlation(high_history, high_current)
    f_low = w1 * horizontal_shift(low_history, low_current, phi) + \
            w2 * scope(low_history, low_current, psi) + \
            w3 * correlation(low_history, low_current)
    f_close = w1 * horizontal_shift(close_history, close_current, phi) + \
              w2 * scope(close_history, close_current, psi) + \
              w3 * correlation(close_history, close_current)
    f_volume = w1 * horizontal_shift(volume_history, volume_current, phi) + \
               w2 * scope(volume_history, volume_current, psi) + \
               w3 * correlation(volume_history, volume_current)
    #
    similarity = min(f_open, f_high, f_low, f_close, f_volume)
    # similarity = (f_open + f_high + f_low + f_close + f_volume)/5
    return similarity


# print(similarity_preference(k_line_full_df[:1], k_line_full_df[:1], 0.2, 0.7, 0.1, 0.5, 0.5))


"""
Looking for the most similar k line
"""


def kline_similar(history_df, current_df, history_stand, current_stand, top_n, method):
    """ looking for the most similar k line(ohlc) in history """
    global distance_i, distance_top_n
    if method == 'DTW':
        # DTW distance
        distance_i = [distance_dtw(history_stand[j:j + 1], current_stand)
                      for j in range(len(history_stand))]
        # print(distance_i)
        distance_top_n = sorted(distance_i)[:top_n]  # top n best distance
    elif method == 'ED':
        distance_i = [ED(history_stand[j:j + 1], current_stand)
                      for j in range(len(history_stand))]
        distance_top_n = sorted(distance_i)[:top_n]  # top n best distance
    elif method == 'preference':
        distance_i = [similarity_preference(history_stand[j:j + 1], current_stand)
                      for j in range(len(history_stand))]
        distance_top_n = sorted(distance_i, reverse=True)[:top_n]  # top n best distance
    # looking for the most similar k line through DTW distance
    ind_top = [distance_i.index(d) for d in distance_top_n]
    k_line_history_top = [history_df[ind:ind + 1] for ind in ind_top]  # top n best history k line
    # get volume data for bar figure
    volume_history_top = [list(kline.volume)[0] for kline in k_line_history_top]
    volume_current_ = list(current_df.volume)[0]
    return k_line_history_top, current_df, volume_history_top, volume_current_, distance_top_n


# s = kline_similar(k_line_full_df, k_line_full_df[1000:1001], k_line_full_df, k_line_full_df[1000:1001],2, 'preference')
# print(s)
# print(k_line_full_df)
# print('hist',s[0])
# print('current',s[1])


"""
Get results
"""


def input(data_stock, L):
    # generate k line data-frame
    date_datetime = [datetime.strptime(d, '%Y/%m/%d') for d in data_stock.date]
    data_stock['date_candlestick'] = [date2num(d) for d in date_datetime]
    data_seg = data_stock.apply(lambda x: seg(x, L), axis=0)  # get time series whose length is L
    data_seg.date = seg(date_datetime, L)
    data_seg = DataFrame({'date': data_seg.date,
                          'date_candlestick': data_seg.date_candlestick,
                          'open': data_seg.open,
                          'high': data_seg.high,
                          'low': data_seg.low,
                          'close': data_seg.close,
                          'volume': data_seg.volume})
    stock = DataFrame({'date': data_seg.date.apply(lambda x: list(x)),
                       'date_candlestick': data_seg.date_candlestick.apply(lambda x: list(x)),
                       'open': data_seg.open.apply(lambda x: list(x)),
                       'high': data_seg.high.apply(lambda x: list(x)),
                       'low': data_seg.low.apply(lambda x: list(x)),
                       'close': data_seg.close.apply(lambda x: list(x)),
                       'volume': data_seg.volume.apply(lambda x: list(x))})  # time series data-frame
    return stock


# data_stock: origin stock data including date(timestamp), ohlc and volume
# parameter L: the length of time series
# start date: the start date of test time series (eg:'2016-4-23')
# top_n: top n best time series
def results(data_stock, start_date, L, top_n, method):
    """ Input: start date, origin stock data, and the length of time series """
    # data standardization
    stock_ohlcv = data_stock[['open', 'high', 'low', 'close', 'volume']]
    data_stock_stand = (stock_ohlcv - mean(stock_ohlcv)) / np.std(stock_ohlcv, ddof=1)
    data_stock_stand['date'] = data_stock['date']

    stock_stand = input(data_stock_stand, L)
    stock = input(data_stock, L)

    # test stock data which start date is 'start date' and history stock data which before the 'start date'
    date_datetime = [datetime.strptime(d, '%Y/%m/%d') for d in data_stock.date]
    date_start = datetime.strptime(start_date, '%Y-%m-%d')  # translate the date into datetime form
    ind_drop = date_datetime.index(date_start)  # the start drop index
    date_drop = date_datetime[ind_drop: ind_drop + L]  # drop indexes and we compare the kline in history data
    # print(date_stamp2time)
    v_bool_test = stock.apply(lambda x: x.date == date_drop, axis=1)
    # print(v_bool_test)
    ind_drop_start = v_bool_test[v_bool_test == True].index.tolist()[0]

    # input data for kline similar function
    k_line_test = stock[v_bool_test]
    k_line_test_stand = stock_stand[v_bool_test]
    k_line_history = stock.drop(range(ind_drop_start - L + 1, len(stock)))
    # print('hist',k_line_history.date[len(k_line_history)-1:len(k_line_history)][0])
    k_line_history_stand = stock_stand.drop(range(ind_drop_start - L + 1, len(stock_stand)))
    res = kline_similar(k_line_history, k_line_test, k_line_history_stand,
                        k_line_test_stand, top_n, method)
    return res


data_s = data_origin_[['date', 'open', 'high', 'low', 'close', 'volume']]
# print(data_s.date)
res = results(data_s, '2017-3-22', 75, 3, 'ED')
print('hist', list(res[0][0].date))
print(list(res[0][1].date))
print(list(res[0][2].date))
print('distance', res[4])

"""
Plot k line figure by candlestick_ohlc
"""


def fig_candlestick(data_stock, start_date, L, top_n, method):
    """ plot k line comparision figure """
    # get similarity results
    res = results(data_stock, start_date, L, top_n, method)
    k_line_history_top = res[0]
    k_line_current = res[1]
    volume_history_top = res[2]
    volume_current = res[3]
    cor = res[4]
    # print(cor)

    # print([(i,j,k,l,m) for i,j,k,l,m in k_line_history_candle_top[0]])
    k_line_current_candle = zip(reduce(merge, k_line_current.date_candlestick),
                                reduce(merge, k_line_current.open),
                                reduce(merge, k_line_current.high),
                                reduce(merge, k_line_current.low),
                                reduce(merge, k_line_current.close))
    # print([(i,j,k,l,m) for i,j,k,l,m in k_line_current_candle])

    # plot k line figures and volume figures
    plt.figure()  # generate figure 1 for test data
    ax1 = plt.subplot(211)  # sub figure1 in figure 1 for ohlc
    ax1.xaxis_date()
    plt.xticks()  # x ticks rotate 45 degree
    plt.yticks()
    plt.title('From' + start_date + 'begin')
    plt.xlabel("Date")
    plt.ylabel("Price")
    candlestick_ochl(ax1, k_line_current_candle, width=0.3, colorup='r', colordown='g')
    plt.subplot(212)  # sub figure2 in figure1 for volume
    plt.bar(range(len(volume_current)), volume_current, 0.4, color='b')

    for i in range(top_n):
        kline = k_line_history_top[i]
        k_line_history_candle_i = zip(reduce(merge, kline.date_candlestick),
                                     reduce(merge, kline.open),
                                     reduce(merge, kline.high),
                                     reduce(merge, kline.low),
                                     reduce(merge, kline.close))
        # print([(i, j, k, l, m) for i, j, k, l, m in k_line_history_candle_i])
        plt.figure()  # generate figure2 for top n best history data
        ax2 = plt.subplot(211)  # sub figure1 in figure 2 for ohlc
        ax2.xaxis_date()
        plt.xticks()
        plt.yticks()
        plt.title('Back from' + start_date)
        plt.xlabel("Date")
        plt.ylabel("Price")
        candlestick_ochl(ax2, k_line_history_candle_i, width=0.3, colorup='r', colordown='g')
        plt.subplot(212)
        plt.bar(range(len(volume_history_top[i])), volume_history_top[i], 0.4, color='b')
        plt.show()


# fig_candlestick(data_s, '2016-4-21', 60, 3, 'preference')
# print(data_origin)


