import processes as proc
import molecules as mol
import replication05 as rep

class Model(object):
    """
    Initializes the states and processes for the model and lets the processes update their corresponding states.
    """
    def __init__(self):

        DNA_length=10000
        ATP_mol=600000
        NT_mol=1400000
        self.states = {}
        self.processes = {}

        # initiate states
        self.Helicase = {'Helicase_{0}'.format(i): rep.Helicase(i, 'Helicase_{0}'.format(i)) for i in xrange(0,2)}
        self.PolymeraseIII = {'PolymeraseIII_{0}'.format(i):rep.PolymeraseIII(i,'PolymeraseIII_{0}'.format(i)) for i in xrange(0,4)}
        self.DNA = {'DNA':rep.DNA('DNA','DNA', DNA_length, 2*DNA_length)}
        self.states.update(self.Helicase)
        self.states.update(self.PolymeraseIII)
        self.states.update(self.DNA)

        # initiate processes
        Replication = rep.Replication(1, "Replication", ATP_mol, NT_mol)
        Replication.set_states(self.DNA.keys(), self.PolymeraseIII.keys() + self.Helicase.keys())
        self.processes = {"Replication":Replication}

    def step(self):
        """
        Do one update step for each process.

        """
        for p in self.processes:
            self.processes[p].update(self)

    def simulate(self, steps, log=True):
        """
        Simulate the model for some time.

        """
        for s in xrange(steps):
            self.step()
            #print (self.Helicase['Helicase_0'].position, self.Helicase['Helicase_0'].bound, self.PolymeraseIII['PolymeraseIII_0'].position, self.PolymeraseIII['PolymeraseIII_0'].bound,self.Helicase['Helicase_1'].position, self.Helicase['Helicase_1'].bound, self.PolymeraseIII['PolymeraseIII_1'].position, self.PolymeraseIII['PolymeraseIII_1'].bound) 

if __name__ == "__main__":
    c = Model()
    c.simulate(150, log=True)
    
