import numpy

def DGEMM_symul(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc):

    if not isinstance(a,int) and len(a.shape)>2:
        a=numpy.squeeze(a)
    if not isinstance(b,int) and len(b.shape)>2:
        b=numpy.squeeze(b)
    if not isinstance(c,int) and len(c.shape)>2:
        c=numpy.squeeze(c)
    #print "A", a.shape
    #print a
    #print "B", b.shape
    #print b
    if transa.upper()=='N':
        if transb.upper()=='N':
            res=alpha*numpy.dot(a[0:m,0:k],b[0:k,0:n])
        else:
            res=alpha*numpy.dot(a[0:m,0:k],b[0:n,0:k].transpose())
    else:
        if transb.upper()=='N':
            res=alpha*numpy.dot(a[0:k,0:m].transpose(),b[0:k,0:n])
        else:
            res=alpha*numpy.dot(a[0:k,0:m].transpose(),b[0:n,0:k].transpose())

    if beta!=0:
        res=res+beta*c

    return res

        
def countCov(x,ip,iwhat):

    sdt=x.shape
    m=sdt[1]
    n=sdt[-1]
    if len(sdt)>2:
        trials=sdt[0]
    else:
        trials=1

    mip=m*ip
    nnconst=n-2*ip

    rleft=numpy.zeros((mip,mip))
    rright=numpy.zeros((mip,m))
    r=numpy.zeros((m,m))
    rleft_tot=numpy.zeros((mip,mip))
    rright_tot=numpy.zeros((mip,m))
    r_tot=numpy.zeros((m,m))

    for trial in range(trials):
        rconst=numpy.zeros((ip+1,m,m))
        for k in range(ip+1):
            #print 0,m,'-',ip-k,ip-k+nnconst,'--',0,m,'-',ip,ip+nnconst
            rconst[k,:,:]=DGEMM_symul('N','T',m,m,nnconst,1,x[trial,0:m,ip-k:ip-k+nnconst],m,x[trial,0:m,ip:ip+nnconst],m,0,0,m)
            if iwhat==8:
                #print 0,m,'-',ip+k,ip+k+nnconst,'--',0,m,'-',ip,ip+nnconst
                rconst[k,:,:]=DGEMM_symul('N','T',m,m,nnconst,1,x[trial,0:m,ip+k:ip+k+nnconst],m,x[trial,0:m,ip:ip+nnconst],m,1,rconst[k,:,:],m)

        for i in range(ip+1):
            for j in range(i,ip+1):
                k=j-i
                r[:,:]=rconst[k,:,:].squeeze()
                if i>0:
                    nn_beg=i
                    #print 0,m,'-',ip-i-k,ip-i-k+nn_beg,'--',0,m,'-',ip-i,ip-i+nn_beg
                    r[:,:]=DGEMM_symul('N','T',m,m,nn_beg,1,x[trial,0:m,ip-i-k:ip-i-k+nn_beg],m,x[trial,0:m,ip-i:ip-i+nn_beg],m,1,r,m)
                if i!=ip:
                    nn_end=ip-i
                    #print 0,m,'-',n-ip-k,n-ip-k+nn_end,'--',0,m,'-',n-ip,n-ip+nn_end
                    r[:,:]=DGEMM_symul('N','T',m,m,nn_end,1,x[trial,0:m,n-ip-k:n-ip-k+nn_end],m,x[trial,0:m,n-ip:n-ip+nn_end],m,1,r,m)

                if iwhat==8:
                    if i<ip:
                        nn_beg=ip-i
                        #print 0,m,'-',i+k,i+k+nn_beg,'--',0,m,'-',i,i+nn_beg
                        r[:,:]=DGEMM_symul('N','T',m,m,nn_beg,1,x[trial,0:m,i+k:i+k+nn_beg],m,x[trial,0:m,i:i+nn_beg],m,1,r,m)
                    if i!=0:
                        nn_end=i
                        #print 0,m,'-',n-ip+k,n-ip+k+nn_end,'--',0,m,'-',n-ip,n-ip+nn_end
                        r[:,:]=DGEMM_symul('N','T',m,m,nn_end,1,x[trial,0:m,n-ip+k:n-ip+k+nn_end],m,x[trial,0:m,n-ip:n-ip+nn_end],m,1,r,m)

                if (i==0) and (j>0):
                    rright[(j-1)*m:j*m,:]=r[:,:]
                if (i>0) and (j>0):
                    rleft[(j-1)*m:(j)*m,(i-1)*m:(i)*m]=r[:,:]
                    rleft[(i-1)*m:(i)*m,(j-1)*m:(j)*m]=r[:,:]

        for column in range(m):
            r[:,column]=rconst[0,:,column]

        nn_end=ip
        #print 0,m,'-',n-ip,n-ip+nn_end,'--',0,m,'-',n-ip,n-ip+nn_end
        r[:,:]=DGEMM_symul('N','T',m,m,nn_end,1,x[trial,0:m,n-ip:n-ip+nn_end],m,x[trial,0:m,n-ip:n-ip+nn_end],m,1,r,m)
        if iwhat==8:
            nn_beg=ip
            #print 0,m,'-',0,nn_beg,'--',0,m,'-',0,nn_beg
            r[:,:]=DGEMM_symul('N','T',m,m,nn_beg,1,x[trial,0:m,0:nn_beg],m,x[trial,0:m,0:nn_beg],m,1,r,m)

        rleft_tot=rleft_tot+rleft
        rright_tot=rright_tot+rright
        r_tot=r_tot+r

    if trials>1:
        rleft_tot=rleft_tot/trials
        rright_tot=rright_tot/trials
        r_tot=r_tot/trials

    if iwhat==4:
        covscale=1./(n-ip)
    elif iwhat==8:
        covscale=1./(2*(n-ip))

    rleft_tot=covscale*rleft_tot
    rright_tot=covscale*rright_tot
    r_tot=covscale*r_tot

    return (rleft_tot,rright_tot,r_tot)


def countCorr(x,ip,iwhat):

    sdt=x.shape
    m=sdt[1]
    n=sdt[-1]
    if len(sdt)>2:
        trials=sdt[0]
    else:
        trials=1

    mip=m*ip

    rleft=numpy.zeros((mip,mip))
    rright=numpy.zeros((mip,m))
    r=numpy.zeros((m,m))
    rleft_tot=numpy.zeros((mip,mip))
    rright_tot=numpy.zeros((mip,m))
    r_tot=numpy.zeros((m,m))

    for trial in range(trials):
        for k in range(1,ip+1):
            if iwhat==1:
                corrscale=1./n
                nn=n-k
                #print 0,m,'-',0,nn,'--',0,m,'-',k,k+nn
                #print 1,m,nn,k,k+nn
                r[:,:]=DGEMM_symul('N','T',m,m,nn,corrscale,x[trial,:,0:nn],m,x[trial,:,k:],m,0,0,m)
                #print r
            elif iwhat==2:
                corrscale=1./(n-k)
                nn=n-k
                #print 0,m,'-',0,nn,'--',0,m,'-',k,k+nn
                r[:,:]=DGEMM_symul('N','T',m,m,nn,corrscale,x[trial,0:m,0:nn],m,x[trial,0:m,k:k+nn],m,0,0,m)

            #print 2,k,m,(k-1)*m,k*m  
            rright[(k-1)*m:k*m,:]=r[:,:]

            if k<ip:
                for i in range(ip-k):
                    #print 3,k,m,(k+i)*m,(k+i+1)*m,i*m,(i+1)*m
                    rleft[(k+i)*m:(k+i+1)*m,i*m:(i+1)*m]=r[:,:]
                    #print 4,k,m,i*m,(i+1)*m,(k+i)*m,(k+i+1)*m
                    rleft[i*m:(i+1)*m,(k+i)*m:(k+i+1)*m]=r[:,:].transpose()


        corrscale=1./n
        #print 5,m,n
        r[:,:]=DGEMM_symul('N','T',m,m,n,corrscale,x[trial,:,:],m,x[trial,:,:],m,0,0,m)

        for k in range(ip):
            #print 6,k,m,k*m,(k+1)*m,k*m,(k+1)*m
            rleft[k*m:(k+1)*m,k*m:(k+1)*m]=r[:,:]

        rleft_tot=rleft_tot+rleft
        rright_tot=rright_tot+rright
        r_tot=r_tot+r

    if trials>1:
        rleft_tot=rleft_tot/trials
        rright_tot=rright_tot/trials
        r_tot=r_tot/trials

    return (rleft_tot,rright_tot,r_tot)


def mult_AR(dat,p,meth_num):

    if len(dat.shape)>2:
        trials=dat.shape[0]
    else:
        trials=1
        dat=numpy.reshape(dat,(1,)+dat.shape)
    chans=dat.shape[1]

    if meth_num<3:
        (rleft,rright,r)=countCorr(dat,p,meth_num)
    else:
        (rleft,rright,r)=countCov(dat,p,meth_num)

    xres=numpy.linalg.solve(rleft,rright)
    mip=chans*p
    Vr=DGEMM_symul('T','N',chans,chans,mip,-1,xres,mip,rright,mip,1,r,chans)

    AR=xres.reshape((p,chans,chans)) # w Pythonie wymiary od 3 dokladaja sie na poczatku
    ARr=AR.transpose((0,2,1))

    return (ARr,Vr)
