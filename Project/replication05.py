import random 
import numpy.random as nr
import processes as proc
import molecules as mol
from processes import Process
from molecules import BioMoleculeCount
from molecules import BioMolecule

'''
Klassen, mit Objekten (DNA, Helicase, Polymerase), Prozess (Replikation). Elternklasse (z.B. BioMolecule, BioMoleculeCount). 
Elternklasse vererbt
'''
class DNA(BioMolecule):
    def __init__(self, mid, name, length, nucleotides, mass=0):
        super(DNA, self).__init__(mid, name, mass)
        self._length=length/2
        self._nucleotides = nucleotides

    @property
    def nucleotides(self):
        return self._nucleotides

    @nucleotides.setter
    def nucleotides(self, value):
        self._nucleotides = value

    @property
    def length(self):
        return self._length

    @length.setter
    def length(self, value):
        self._length = value

class PolymeraseIII(BioMoleculeCount):
	#position=int; bound=boolean
    def __init__(self, mid, name, count=0, position=0, bound=False):
        super(PolymeraseIII, self).__init__(mid, name, count)
        self._position=position
        self._bound=bound

    @property
    def bound(self):
        return self._bound

    @bound.setter
    def bound(self, value):
        self._bound = value

    @property
    def position(self):
        return self._position

    @position.setter
    def position(self, value):
        self._position = value

class Helicase(BioMoleculeCount):
	#position=int; bound=boolean
    def __init__(self, mid, name, count=0, position=0, bound=False):
        super(Helicase, self).__init__(mid, name, count)
        self._position=position
        self._bound = bound

    @property
    def bound(self):
        return self._bound

    @bound.setter
    def bound(self, value):
        self._bound = value

    @property
    def position(self):
        return self._position

    @position.setter
    def position(self, value):
        self._position = value


class Replication(Process):

    def __init__(self, id, name, ATP, NT, double=False):
        super(Replication, self).__init__(id, name)
        self._double = double
        self.ATP_molecules = ATP
        self.Nucleotide = NT
        PolymeraseIII_clock_lag_double = 0
        PolymeraseIII_anticlock_lag_double =0

    @property
    def double(self):
        return self._double

    @double.setter
    def double(self, value):
        self._double = value

    

    def update(self, model):
    	#return position helicase und polymeraseIII
        #Uebergeben der Objekte, vom Ori aus wird in zwei Richtungen repliziert (anticlock und clock)
        self.PolymeraseIII_anticlock_lead = model.states['PolymeraseIII_0']
        self.PolymeraseIII_anticlock_lag = model.states['PolymeraseIII_1']
        self.PolymeraseIII_clock_lead = model.states['PolymeraseIII_2']
        self.PolymeraseIII_clock_lag = model.states['PolymeraseIII_3']
        self.Helicase_anticlock = model.states['Helicase_0']
        self.Helicase_clock = model.states['Helicase_1']
        self.DNA = model.states['DNA']


        if self.double ==False: # Replikation darf nur stattfinden, wenn DNA noch nicht veerdoppelt ist
        #in anticlock_Richtung
            if (self.PolymeraseIII_anticlock_lead.bound == False) and (self.PolymeraseIII_anticlock_lag.bound == False) \
            and(self.PolymeraseIII_anticlock_lead.position <self.DNA.length) and(self.PolymeraseIII_anticlock_lag.position <self.DNA.length):
                self.Helicase_anticlock=self.initiate(DNA, self.Helicase_anticlock, self.PolymeraseIII_anticlock_lead, self.PolymeraseIII_anticlock_lag)
            elif (self.Helicase_anticlock.bound == True) and (self.Helicase_anticlock.position >=1500):
                if self.PolymeraseIII_anticlock_lead.bound == True:
                    self.Helicase_anticlock, self.PolymeraseIII_anticlock_lead=self.elongate_lead(DNA, self.PolymeraseIII_anticlock_lead, self.PolymeraseIII_anticlock_lag, self.Helicase_anticlock)   
                if self.PolymeraseIII_anticlock_lag.bound == True:
                    self.PolymeraseIII_anticlock_lag=self.elongate_lag(DNA, self.PolymeraseIII_anticlock_lead, self.PolymeraseIII_anticlock_lag, self.Helicase_anticlock, PolymeraseIII_anticlock_lag_double)
            elif self.Helicase_anticlock.position >=1500:
                self.PolymeraseIII_anticlock_lead=self.terminate_lead(DNA, self.PolymeraseIII_anticlock_lead, self.Helicase_anticlock)
                self.PolymeraseIII_anticlock_lag=self.terminate_lag(DNA, self.PolymeraseIII_anticlock_lag, self.Helicase_anticlock, PolymeraseIII_anticlock_lag_double)    
            if (self.PolymeraseIII_anticlock_lead.bound == False) and (self.PolymeraseIII_anticlock_lead.position <self.DNA.length) and (self.Helicase_anticlock.position >= 1500):
                self.PolymeraseIII_anticlock_lead = self.binding(self.PolymeraseIII_anticlock_lead)
            if (self.PolymeraseIII_anticlock_lag.bound == False) and (self.PolymeraseIII_anticlock_lag.position <self.DNA.length):
                if (self.Helicase_anticlock.position >= 2000):
                    self.PolymeraseIII_anticlock_lag = self.binding(self.PolymeraseIII_anticlock_lag)
                    if self.PolymeraseIII_anticlock_lag.bound == True:
                        
                        if self.Helicase_anticlock.position <3000:
                            self.PolymeraseIII_anticlock_lag.position += 1500
                        else:
                            PolymeraseIII_anticlock_lag_double +=1500 
                            self.PolymeraseIII_anticlock_lag.position += 3000

        #in clock_Richtung
            if (self.PolymeraseIII_clock_lead.bound == False) and (self.PolymeraseIII_clock_lag.bound == False) \
            and(self.PolymeraseIII_clock_lead.position <self.DNA.length) and(self.PolymeraseIII_clock_lag.position <self.DNA.length):
                self.Helicase_clock=self.initiate(DNA, self.Helicase_clock, self.PolymeraseIII_clock_lead, self.PolymeraseIII_clock_lag)
            elif self.Helicase_clock.bound == True :
                if self.PolymeraseIII_clock_lead.bound == True:
                    self.Helicase_clock, self.PolymeraseIII_clock_lead=self.elongate_lead(DNA, self.PolymeraseIII_clock_lead, self.PolymeraseIII_clock_lag, self.Helicase_clock)   
                if self.PolymeraseIII_clock_lag.bound == True:
                    self.PolymeraseIII_clock_lag=self.elongate_lag(DNA, self.PolymeraseIII_clock_lead, self.PolymeraseIII_clock_lag, self.Helicase_clock, PolymeraseIII_clock_lag_double)
            elif self.Helicase_clock.position>=1500:
                self.PolymeraseIII_clock_lead=self.terminate_lead(DNA, self.PolymeraseIII_clock_lead, self.Helicase_clock)
                self.PolymeraseIII_clock_lag=self.terminate_lag(DNA, self.PolymeraseIII_clock_lag, self.Helicase_clock, PolymeraseIII_clock_lag_double)    
            if (self.PolymeraseIII_clock_lead.bound == False) and (self.PolymeraseIII_clock_lead.position <self.DNA.length) and (self.Helicase_clock.position >= 1500):
                self.PolymeraseIII_clock_lead = self.binding(self.PolymeraseIII_clock_lead)
            if (self.PolymeraseIII_clock_lag.bound == False) and (self.PolymeraseIII_clock_lag.position <self.DNA.length) and (self.Helicase_clock.position >= 2000):
                self.PolymeraseIII_clock_lag = self.binding(self.PolymeraseIII_clock_lag)
                if self.PolymeraseIII_clock_lag.bound == True:
                    if self.Helicase_clock.position <3000:
                        self.PolymeraseIII_clock_lag.position += 1500
                    else:
                        PolymeraseIII_clock_lag_double +=1500 
                        self.PolymeraseIII_clock_lag.position += 3000
        #wenn beide Richtungen komplett sind, wird double auf True gesetzt, um weitere Replikation zu verhindern
            if (self.PolymeraseIII_anticlock_lead.position >= self.DNA.length) and (self.PolymeraseIII_anticlock_lag.position >= self.DNA.length) and (self.PolymeraseIII_clock_lead.position >= self.DNA.length) and (self.PolymeraseIII_clock_lag.position >= self.DNA.length):
                self.double = True

        self.DNA.nucleotides = 4*self.DNA.length + self.PolymeraseIII_anticlock_lead.position + self.PolymeraseIII_anticlock_lag.position+ self.PolymeraseIII_clock_lag.position +self.PolymeraseIII_clock_lead.position
        #print self.DNA.nucleotides
        print self.Helicase_anticlock.position, self.PolymeraseIII_anticlock_lead.position, self.PolymeraseIII_anticlock_lag.position,self.Helicase_clock.position, self.PolymeraseIII_clock_lead.position, self.PolymeraseIII_clock_lag.position,self.DNA.nucleotides 
        print self.Helicase_anticlock.bound, self.PolymeraseIII_anticlock_lead.bound, self.PolymeraseIII_anticlock_lag.bound
       

    def initiate(self, DNA,Hel, Pol_lead, Pol_lag):
    	"""
        wird aufgerufen, wenn Polymerase noch nicht gebunden ist, Helikase bindet mit def. Bindungswahrscheinlichkeit und startet 
        Strangauftrennung (und damit ATP-Verbrauch). Falls Abstand Helikase-Polymerase Mindestabstand erreicht hat, bindet Polymerase mit def. Bindungswahrscheinlichkeit.
        """
        Helicase = Hel
        PolymeraseIII_lead =Pol_lead
        PolymeraseIII_lag = Pol_lag
        if Helicase.bound == False:
            x_number =nr.randint(1,10)#Bindungswahrscheinlichkeit Helikase
            if x_number ==1:
                Helicase.bound =True
        elif self.ATP_molecules >= 100 and ((Helicase.position - PolymeraseIII_lead.position) < 3000) and ((Helicase.position - PolymeraseIII_lag.position) < 3000):
            Helicase.position += 100
            self.ATP_molecules -= 100
        elif self.ATP_molecules > 0 and (Helicase.position - PolymeraseIII_lead.position) < 3000 and (Helicase.position - PolymeraseIII_lag.position) < 3000:
            Helicase.position += self.ATP_molecules
            self.ATP_molecules -= self.ATP_molecules
    
        if Helicase.position > self.DNA.length:
            self.ATP_molecules=self.ATP_molecules+(Helicase.position -self.DNA.length)
            Helicase.position = self.DNA.length
        #print ('ATP:',self.ATP_molecules,'NT:',self.Nucleotide)
        return Helicase

    def binding(self, Pol):
        PolymeraseIII=Pol
        y_number =nr.randint(1,5)# Bindungswahrscheinlichkeit Polymerase
        if y_number ==1:
            PolymeraseIII.bound = True
        return PolymeraseIII

    def elongate_lead(self,DNA, Pol_lead, Pol_lag, Hel):
        """
        Wird aufgerufen, wenn Polymerase und Helicase gebunden sind. Testet, ob genug ATP Molekuele und Nukelotide vorhanden sind.
        Verlaengert pro Step um 100 Nukelotide oder die maximal moegliche Anzahl bei ATP/Nukleotid Begrenzung.
        Der maximale Abstand zwischen Helikase und Polymerase ist 3000, der minimale 1500.
        """
        Helicase = Hel
        PolymeraseIII_lead = Pol_lead
        PolymeraseIII_lag = Pol_lag
        if self.ATP_molecules >= 100 and (Helicase.position - PolymeraseIII_lead.position) < 3000 and (Helicase.position - PolymeraseIII_lag.position) < 3000: #genug ATP, Abstand klein genug
            Helicase.position += 100 
            self.ATP_molecules -= 100
            if self.Nucleotide >= 200 and (Helicase.position - PolymeraseIII_lead.position) > 1500: #genug  Nucleotide (>=200)
                PolymeraseIII_lead.position += 100
                self.Nucleotide -= 200
            elif self.Nucleotide > 1 and (Helicase.position - PolymeraseIII_lead.position) > 1500: #nicht genug Nucleotide (1-199)
                PolymeraseIII.position_lead += self.Nucleotide/2
                Helicase.position = Helicase.position -100 +self.Nucleotide/2
                self.ATP_molecules =self.ATP_molecules+100-self.Nucleotide/2
                self.Nucleotide -= 2*(self.Nucleotide/2)
        
        elif self.ATP_molecules >= 0 and (Helicase.position - PolymeraseIII_lead.position) < 3000 and (Helicase.position - PolymeraseIII_lag.position) < 3000: #nicht genug ATP, Abstand klein genug
            Helicase.position += self.ATP_molecules
            if self.Nucleotide >= 200 and (Helicase.position - PolymeraseIII_lead.position) > 1500:  #genug Nucleotide
                PolymeraseIII_lead.position += 100
                self.Nucleotide -= 200
            elif self.Nucleotide > 1 and (Helicase.position - PolymeraseIII_lead.position) > 1500: #nicht genug Nucleotide
                PolymeraseIII_lead.position += self.Nucleotide/2
                Helicase.position = Helicase.position -self.ATP_molecules +self.Nucleotide/2
                self.ATP_molecules -=self.Nucleotide/2
                self.Nucleotide -= 2*(self.Nucleotide/2)
            self.ATP_molecules -= self.ATP_molecules

        if Helicase.position > self.DNA.length:
            self.ATP_molecules=self.ATP_molecules+(Helicase.position -self.DNA.length)
            Helicase.position = self.DNA.length

        if Helicase.position >= self.DNA.length:
            Helicase.bound =False
        #print ('ATP:',self.ATP_molecules,'NT:',self.Nucleotide)
        return Helicase, PolymeraseIII_lead

    def elongate_lag(self, DNA, Pol_lead, Pol_lag, Hel, Pol_lag_double):
        Helicase = Hel
        PolymeraseIII_lead = Pol_lead
        PolymeraseIII_lag = Pol_lag
        PolymeraseIII_lag_double = Pol_lag_double
        
        if self.Nucleotide >= 200: #and (Helicase.position - PolymeraseIII_lag.position) > 2000: #genug  Nucleotide (>=200)
            PolymeraseIII_lag.position -= 100
            self.Nucleotide -= 200
        elif self.Nucleotide > 1: #and (Helicase.position - PolymeraseIII_lag.position) > 2000: #nicht genug Nucleotide (1-199)
            PolymeraseIII_lag.position -= self.Nucleotide/2
            self.Nucleotide -= 2*(self.Nucleotide/2)

        if PolymeraseIII_lag.position <= PolymeraseIII_lag_double:
            PolymeraseIII_lag.bound =False
        #print ('ATP:',self.ATP_molecules,'NT:',self.Nucleotide)
        
        return PolymeraseIII_lag



    def terminate_lead(self,DNA, Pol, Hel):
        """
        Beendet die Replikation. Wird aufgerufen, wenn die Helicase nicht mehr gebunden ist.
        Wenn genug Nucleotide vorhanden sind, wird die Polymerase pro Step um 100 Nukelotide verschoben,
        sonst um die  maximal moegliche Anzahl. Ist die Mitte der DNA erreicht, wird die Polymerase abgeloest
        """
        Helicase = Hel
        PolymeraseIII = Pol
        
        #print self.DNA.length, PolymeraseIII.position
    	#wenn pos= length helicase + polyIII fallen ab
        if PolymeraseIII.bound == True:
            if self.Nucleotide >= 200 :
                PolymeraseIII.position += 100
                self.Nucleotide -= 200
            elif self.Nucleotide>1 :
                PolymeraseIII.position += self.Nucleotide/2
                self.Nucleotide -= 2*self.Nucleotide/2

        if PolymeraseIII.position >= self.DNA.length:
            self.Nucleotide += (PolymeraseIII.position-self.DNA.length)*2
            PolymeraseIII.position=self.DNA.length
            PolymeraseIII.bound = False

        return PolymeraseIII

    def terminate_lag(self,DNA, Pol, Hel, Pol_lag_double):
        """
        Beendet die Replikation. Wird aufgerufen, wenn die Helicase nicht mehr gebunden ist.
        Wenn genug Nucleotide vorhanden sind, wird die Polymerase pro Step um 100 Nukelotide verschoben,
        sonst um die  maximal moegliche Anzahl. Ist die Mitte der DNA erreicht, wird die Polymerase abgeloest
        """
        Helicase = Hel
        PolymeraseIII = Pol
        PolymeraseIII_lag_double = Pol_lag_double
        
        #print self.DNA.length, PolymeraseIII.position
        #wenn pos= length helicase + polyIII fallen ab
        if PolymeraseIII.bound == True:
            if self.Nucleotide >= 200 : #genug  Nucleotide (>=200)
                PolymeraseIII.position -= 100
                self.Nucleotide -= 200
            elif self.Nucleotide > 1 : #nicht genug Nucleotide (1-199)
                PolymeraseIII.position -= self.Nucleotide/2
                self.Nucleotide -= 2*(self.Nucleotide/2)

        if PolymeraseIII.position <= PolymeraseIII_lag_double:
            PolymeraseIII.bound =False

        if PolymeraseIII.position >= self.DNA.length:
            self.Nucleotide += (PolymeraseIII.position-self.DNA.length)*2
            PolymeraseIII.position=self.DNA.length
            PolymeraseIII.bound = False

        return PolymeraseIII

    def gene_check(self,DNA,Pol_ac,Pol_c,gene_begin,gene_end): 
        """
        Testet, ob ein Gen schon doppelt oder einfach vorliegt. Uebergeben werden muessen die DNA-Polymerasen,
        sowie Start und Endpunkt des Gens auf dem Strang.
        return=2 fuer Gene, die schon repliziert wurden
        return=1 fuer noch nicht repliziert wurden
        """
        PolymeraseIII_ac = Pol_ac
        PolymeraseIII_c = Pol_c
        if (gene_end < PolymeraseIII_c.position) or (gene_begin > (2*self.DNA.length-PolymeraseIII_ac.position)):
            return 2
        else:
            return 1
