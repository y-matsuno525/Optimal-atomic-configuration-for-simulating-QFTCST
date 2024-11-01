import math
import numpy as np
from scipy.optimize import root_scalar
import matplotlib.pyplot as plt

#角度はすべてラジアン表記
theta_m=0.9546951 #双極子相互作用が消える方向を表す角度
theta_list_NNN=[] #NNのinternuclear axisの方向を表す角度(NNNのinternuclear axisを基準とした)を格納するリスト
theta_list_q=[] #NNのinternuclear axisの方向を表す角度(NNNのquantization axisを基準とした)を格納するリスト
d_NN_list=[]
d_NNN_list=[]
kappa_list=[] #結合定数のリスト

#初期設定
l=1 #NNNのinternuclear axis間の距離（あとで可変にするかも）
L=601 #サイト数
c=100
d=0.1

#metricがencodeされた関数κの定義
def kappa(n,alpha=10):
    return alpha*np.tanh(((n-(L-1)/2)-0.5)*d)/(4*d)

#配置を決める方程式の定義
def equation(x,n):
    return 2*c*(1-3*np.cos(theta_m-x)**2)+kappa(n)*(l/np.sin(x))**3

# 解を探す関数
def find_root_in_range(func, n,lower, upper):
    try:
        result = root_scalar(func,args=(n,) ,method='brentq',bracket=[lower, upper])
        if result.converged:
            return result.root
        else:
            return None
    except ValueError:
        return None

for n in range(L):
    if n==0:
        root = find_root_in_range(equation, n,0.0001, np.pi-theta_m)
        if root is not None:
            theta_list_q.append(theta_m+root) 
            theta_list_NNN.append(root)
            d_NN_list.append(l/math.sin(theta_list_NNN[0]))
            

        else:
            print(f"Root not found for n={n} in initial case.")
            break
    else:
        if n%2 != 0:
            lower_bound=0.0001
            upper_bound=(np.pi+theta_m-theta_list_q[n-1])/2
            if upper_bound > theta_m: #NNNとNNの関係を維持するための条件
                upper_bound = theta_m
            root = find_root_in_range(equation, n,lower_bound, upper_bound) #解が複数ある場合は？
            if root is not None:
                theta_list_NNN.append(root)
                theta_list_q.append(theta_m-root)
                d_NN_list.append(l/math.sin(theta_list_NNN[n]))
            else:
                print(f"Root not found for n={n} in even case.")
                break
        else:
            lower_bound=0.0001
            upper_bound=np.pi-2*theta_m+2*theta_list_q[n-1]
            if upper_bound > theta_m:
                upper_bound = theta_m
            root = find_root_in_range(equation, n,lower_bound, upper_bound) #解が複数ある場合は？
            if root is not None:
                theta_list_NNN.append(root)
                theta_list_q.append(root+theta_m)
                d_NN_list.append(l/math.sin(theta_list_NNN[n]))
            else:
                print(f"Root not found for n={n} in odd case.")
                break
    kappa_list.append(-2*c*(1-np.cos(theta_list_q[n])**2)/d_NN_list[n]**3)

#a=[x*180/np.pi for x in theta_list_q]
#print(a)
#print(d_NN_list)

#plt.figure(figsize=(8,7))
#plt.plot(kappa_list)
#plt.show()