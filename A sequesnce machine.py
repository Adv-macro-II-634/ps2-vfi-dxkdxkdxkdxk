A = [1.1,0.678]

transmatrix = [[0.977,1-0.977],[1-0.926,0.926]]

def forecast(days): 
    
    Today = '1.1'
    #print("start state" + Today)
    actlist = [Today]
    i = 0
    #prob = 1 
    
    while i != days:
        if Today == '0.678':
          change = np.random.choice(A,p=transmatrix[1])
        
          if change == 0.678:
          # prob = prob * 0.926
            actlist.append('0.678')
            pass
        
          else:
           # prob = prob * (1-0.926)
            actlist.append('1.1')
            
   
        else:
          change = np.random.choice(A,p=transmatrix[0])
          if change == 1.1:
          #  prob = prob * 0.977
            actlist.append('1.1')
            pass
        
          else:
            #prob = prob * (1-0.977)
            actlist.append('0.678')
            
        i += 1
    return actlist

forecast(1000)
