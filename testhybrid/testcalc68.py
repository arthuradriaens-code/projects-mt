import numpy as np

#hist = [0.025,0.1,0.65, 0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.025,0.025,0.025,0.025,0.075,0.05 ]
#bins = [-0.32437862,-0.1696811,-0.01498358 ,0.13971394 ,0.29441146 ,0.44910898
#,0.6038065 ,0.75850401,0.91320153,1.06789905,1.22259657,1.37729409
#,1.53199161,1.68668913,1.84138665,1.99608417,2.15078169,2.3054792
#,2.46017672,2.61487424,2.76957176]
hist = np.array([0,0,0,0.69,0,0.31])
bins = [-2,-1,0,0.099,0.1,0.2,0.3]

def Calc68(hist,bins):
    totalleft = 0
    totalright = 0
    for i,histpoint in enumerate(hist):
        totalleft += histpoint
        if totalleft >= 0.16:
            leftindex = i
            break
    for i,histpoint in enumerate(hist[::-1]):
        totalright += histpoint
        if totalright >= 0.16:
            rightindex = i
            break
    rightindex = len(hist)-rightindex-1
    print(hist[leftindex:rightindex].sum())
    while hist[leftindex:rightindex].sum() >= 0.68:
        leftindex += 1
    leftindex -=1
    while hist[leftindex:rightindex].sum() >= 0.68:
        rightindex -= 1
    rightindex +=1
    print(hist[leftindex:rightindex].sum())
    print(leftindex)
    print(rightindex)
    print("interval =" + str(bins[rightindex] - bins[leftindex]))


        
Calc68(hist,bins)
