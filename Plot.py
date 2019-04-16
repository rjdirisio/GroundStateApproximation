import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

class myPlot:
    def __init__(self,cfg,typ,rangey,rangex,xl,yl,overlay,indv,pltdata,binz1,binz2=[],fancy=False):
        self.cfg = cfg
        self.type = typ
        self.rangey = rangey
        self.rangex = rangex
        self.xl = xl
        self.yl = yl
        self.overlay = overlay
        self.indv = indv
        self.pltdata = pltdata
        self.binCenters1 = binz1 #binz no longer is just a number, but an array of bin centers
        self.binCenters2 = binz2  # binz no longer is just a number, but an array of bin centers
        self.fancy = fancy
    def determineSmallestMax2D(self):
        overall = 1000000000.
        for i in self.pltdata:
            thisMax = np.amax(i)
            if thisMax < overall:
                overall = thisMax
        return overall

    def getYview(self):
        overall = 0.
        for i in self.pltdata:
            thisMax = np.amax(i)
            if thisMax > overall:
                overall = thisMax
        return overall

    def plotIt(self):
        titlemod = ''
        if self.type == '1d':
            if self.indv == True: #Standard 1dHists
                print self.xl,self.yl
                plt.xlabel(self.xl)
                plt.ylabel(self.yl)
                plt.plot(self.binCenters1,self.pltdata,'k',linewidth=2)
                plt.ylim(bottom=0.0)
                #plt.show()
                savef = titlemod+self.yl.replace(" ", "")+'vs'+self.xl.replace(" ", "")+self.cfg+'.png'
                plt.savefig(savef)
            elif self.indv == False:
                num=1
                if self.overlay == True: #self.overlay 1dHists
                    for i in self.pltdata:
                        plt.plot(self.binCenters1,i,lienwidth=2,label=num)
                        num+=1
                    plt.xlabel(self.xl)
                    plt.ylabel(self.yl)
                    plt.ylim(bottom=0.0)
                    plt.legend()
                    savef = titlemod+self.yl.replace(" ", "")+'vs'+self.xl.replace(" ", "")+self.cfg+'Overlay.png'
                    plt.savefig(savef)
                    
                elif self.overlay == False: #self.indv Files
                    num=0
                    yv = self.getYview()
                    for i in self.pltdata:
                        plt.figure()
                        plt.plot(self.binCenters1,i)
                        plt.xlabel(self.xl)
                        plt.ylabel(self.yl)
                        plt.ylim(0,yv)
                        plt.title("WF "+str(num+1))
                        savef = titlemod+self.yl.replace(" ", "")+'vs'+self.xl.replace(" ", "")+str(num)+self.cfg+'.png'
                        num+=1
                        plt.savefig(savef)
                        
        elif self.type == '2d':
            if self.indv == True: #Standard 1dHists
                if not self.fancy:
                    plt.contour(self.binCenters1,self.binCenters2,self.pltdata,colors='k')
                    plt.contourf(self.binCenters1,self.binCenters2,self.pltdata)
                    plt.colorbar()

                else:
                    print 'white BG activated)'
                    cmap = mpl.cm.Greys
                    norm = mpl.colors.Normalize(vmin=np.amin(self.pltdata) + 3, vmax=maxDat)
                    ax = plt.contourf(self.binCenters1, self.binCenters2, i, norm=norm, cmap=cmap)
                    ax.cmap.set_under("w")
                    plt.colorbar(ax)
                plt.xlabel(self.xl)
                plt.ylabel(self.yl)
                savef = self.yl.replace(" ", "")+'vs'+self.xl.replace(" ", "")+self.cfg+'.png'
                plt.savefig(savef)

            elif self.indv==False:
                if self.overlay==False:
                    maxDat = self.determineSmallestMax2D()
                    num=1
                    print np.size(self.pltdata)
                    print np.size(self.pltdata[0])
                    print np.size(self.pltdata[:][:][0])
                    
                    for i in self.pltdata:
                        self.pltdata[self.pltdata>maxDat] = maxDat
                        plt.contour(self.binCenters1,self.binCenters2,i,colors='k')
                        plt.contourf(self.binCenters1,self.binCenters2,i)
                        plt.colorbar()
                        plt.title('WF '+str(num))
                        plt.xlabel(self.xl)
                        plt.ylabel(self.yl)
                        plt.legend()
                        savef = self.yl.replace(" ", "")+'vs'+self.xl.replace(" ", "")+str(num)+self.cfg+'.png'
                        num+=1
                        plt.savefig(savef)
                        plt.close()
        plt.close()


#Examples of how to Plot Stuff
    #PltHists1D(thing,bound,title,xl,yl,overly):
    #PltHists1D(cfg,OuterWaterAngle,(50,150),'Outer HOH Angle','Probability Density',False,weightArray)
    #PltHists1D(cfg,spAngle,(50,150),'Outer HOH Angle', 'Probability Density',False,spweits)
    #PltHists2D(cfg,OuterWaterAngle,InnerWaterAngle,(90,140),(90,140),'Inner HOH Angle','Outer HOH Angle',False,weightArray)
    #PltHists2D(cfg,spOuterAngle,spInnerAngle,(90,140),(90,140),'Inner HOH Angle','Outer HOH Angle',False,spweits)


