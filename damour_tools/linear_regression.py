# Linear regression to find a and b such as :
#  ! y = a*x + b
#
#  ! a = \sum_i (x_i - \bar{x})(y_i - \bar{y}) / \sum_i (x_i - \bar{x})^2
#  ! b = \bar_y - a \bar{x}
#
#  ! \bar{x} = 1/ n \sum_i x_i 
#  ! \bar{y} = 1/ n \sum_i y_i 
#
#  ! R^2 = 1 - SSR / SST
#
#  ! SSR = \sum_i (y_i - \hat{y})^2
#  ! SST = \sum_i (y_i - \bar{y})^2
#
#  ! \hat{y}_i = a*x_i + b


def lin_reg(x,y,nb_points):
    tmp_x=[]
    for i in range(nb_points):
        tmp_x.append(x[len(x)-1-i])

    tmp_y=[]
    for i in range(nb_points):
        tmp_y.append(y[len(y)-1-i])

    bar_x = 0.0
    for i in range(len(tmp_x)):
        bar_x = bar_x + tmp_x[i]
    bar_x = bar_x / len(tmp_x)
 
    bar_y = 0.0
    for i in range(len(tmp_y)):
        bar_y = bar_y + tmp_y[i]
    bar_y = bar_y / len(tmp_y)

    a1 = 0.0
    for i in range(len(tmp_x)):
        a1 = a1 + (tmp_x[i]-bar_x)*(tmp_y[i]-bar_y)
    a2 = 0.0
    for i in range(len(tmp_x)):
        a2 = a2 + (tmp_x[i]-bar_x)**2
    
    a = a1/a2

    b = bar_y - a * bar_x

    hat_y = []
    for i in range(len(tmp_x)):
        hat_y.append(a * tmp_x[i] + b) 
    
    SSR = 0.0
    for i in range(len(tmp_x)):
        SSR = SSR + (tmp_y[i] - hat_y[i])**2

    SST = 0.0
    for i in range(len(tmp_x)):
        SST = SST + (tmp_y[i] - bar_y)**2

    R2 = 1.0 - SSR/SST
 
    return(a,b,R2)

def lin_reg_v2(x,y,weight,nb_points):
    import numpy as np  

    # last n points and reverse their order
    tmp_x = np.array(x[-nb_points:])
    tmp_x = tmp_x[::-1]

    tmp_y = np.array(y[-nb_points:])
    tmp_y = tmp_y[::-1]

    tmp_w = np.array(weight[-nb_points:])
    tmp_w = tmp_w[::-1]

    # linear regression
    fit = np.polynomial.polynomial.polyfit(tmp_x, tmp_y, deg=1, rcond=None, full=False, w=tmp_w)

    # f(x) = ax + b
    a = fit[1] 
    b = fit[0]

    # bar_y = \sum_i^n tmp_y[i]/n
    bar_y = np.sum(tmp_y) / len(tmp_y)

    # hat_y[i] = a*tmp_x[i] + b 
    hat_y = a * tmp_x + b
   
    # SSR = \sum_i (tmp_y[i] - hat_y[i])**2
    SSR = np.sum((tmp_y-hat_y)**2)

    # SST = \sum_i (tmp_y[i] - bar_y)**2
    SST = np.sum((tmp_y - bar_y)**2)

    R2 = 1.0 - SSR/SST

    return(a,b,R2)

# for org mode
#
#+NAME: toto
#+BEGIN_SRC bash :result replace :exports none
#cat file.fci.out.dat
#+END_SRC
#
#+BEGIN_SRC python :var data=toto :result value :exports results
#from linear_regression import create_extrapolation_array
#
#res = create_extrapolation_array(data,1,"PT2")
#return res
#+END_SRC

def create_extrapolation_array(data,state,kind):
    n = len(data)-1

    res = []

    text = ("n","coef","extrapolation (E_h)","R^2")
    res.append(text)

    column = (state-1)*5+1

    for nb_points in range(3,8):
        E = []
        for k in range(n):
            E.append(data[k+1][column])

        PT2 = []
        for k in range(n):
            PT2.append(data[k+1][column+1])

        rPT2 = []
        for k in range(n):
            rPT2.append(data[k+1][column+3])

        if (kind == "PT2"):
            toto=lin_reg(PT2,E,nb_points)
        else:
            toto=lin_reg(rPT2,E,nb_points)

        toto=(nb_points,)+toto
        res.append(toto)

    return(res)
