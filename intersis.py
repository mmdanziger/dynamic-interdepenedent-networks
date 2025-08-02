import scipy.sparse as sps
import networkx as nx
from scipy.integrate import odeint
import numpy as np

def generate_scale_free(N,gamma):
    return nx.barabasi_albert_graph(N,5)
    while True:
        s = []
        while len(s) < N:
            nextval = int(nx.utils.powerlaw_sequence(1, gamma)[0])
            if nextval != 0 and nextval < np.sqrt(N):
                s.append(nextval)
        if sum(s) % 2 == 0:
            break
    G = nx.random_degree_sequence_graph(s)
    G = nx.Graph(G)  # remove parallel edges
    G.remove_edges_from(G.selfloop_edges())
    return G

all_solvers = ["RK45", "RK23", "Radau", "BDF", "LSODA"]
class InterSIS(object):
    def __init__(self,N,k,beta,dist="sf"):
        self.N=N
        self.k=k
        self.beta=beta
        gamma = k
        if dist == "sf":
            self.G1 = generate_scale_free(N,gamma)
            self.G2 = generate_scale_free(N,gamma)
        elif dist == "er":
            self.G1 = nx.erdos_renyi_graph(N,k/N)
            self.G2 = nx.erdos_renyi_graph(N,k/N)
        elif dist == "rr":
            self.G1 = nx.random_regular_graph(k,N)
            self.G2 = nx.random_regular_graph(k, N)
        self.A1 = nx.to_scipy_sparse_matrix(self.G1,format="csr")
        self.A2 = nx.to_scipy_sparse_matrix(self.G2,format="csr")
        self.x1 = np.zeros((N,1))
        self.x2 = np.zeros((N,1))
        self.interaction_type="competitive"
        self.f = np.ones((N,1))
        self.antif = 1 - self.f
        self.hybrid = True
        self.degvec1 = np.array(self.A1.sum(axis=0)).reshape(self.x1.shape)
        self.degvec2 = np.array(self.A2.sum(axis=0)).reshape(self.x2.shape)
        self.degsum1 = sum(self.degvec1.flatten())
        self.degsum2 = sum(self.degvec2.flatten())
        self.drivers = False
        self.diffusion_mode = False
        self.local_mf= True
        self.degree_weighted_op = True
        self.solver = all_solvers[1]

    def set_f(self,f):
        self.f = np.array( np.random.rand(self.N,1) < f).astype(int)
        self.antif = 1  - self.f

    def add_drivers(self,driver_fraction=0.1,source_strength=1,sink_strength=1):
        self.drivers = True
        self.x1driver = np.array( np.random.rand(self.N,1) < driver_fraction).astype(int)
        self.x2driver = np.array(np.random.rand(self.N, 1) < driver_fraction).astype(int)
        self.ddriver1 = source_strength * self.x1driver - sink_strength * (1 - self.x1driver)
        self.ddriver2 = source_strength * self.x2driver - sink_strength * (1 - self.x2driver)

    def apply_to_f(self,x):
        return self.f * x + self.antif

    def interaction(self,x,interaction_type=None):
        if not interaction_type:
            interaction_type = self.interaction_type
        if interaction_type == "interdependent":
            return x
        elif interaction_type == "competitive":
            return 1 - x
        else:
            raise(ValueError)

    def randinit_upto(self,p):
        self.x1 = p[0]*np.random.rand(self.N,1)
        self.x2 = p[1]*np.random.rand(self.N,1)
        if self.drivers:
            for i in range(self.N):
                if self.x1driver[i] == 1:
                    self.x1[i] = 1
                if self.x2driver[i] == 1:
                    self.x2[i] = 1

    def inter_term(self,x,one_or_two):
        if self.local_mf:
            if one_or_two == 1:
                out = self.A1.dot(x) / self.degvec1
            elif one_or_two == 2:
                out = self.A2.dot(x) / self.degvec2
            else:
                raise ValueError
        else:
            out = x
        return np.array(np.nan_to_num(out))


    def order_param(self,x,one_or_two,global_op_override=False):
        if self.degree_weighted_op and not global_op_override:

            if one_or_two == 1:
                out = sum(self.A1.dot(x)) / self.degsum1
            elif one_or_two == 2:
                out = sum(self.A2.dot(x)) / self.degsum2
            else:
                raise ValueError
            return np.array(np.nan_to_num(out)).flatten()
        else:
            return np.mean(x,axis=0)

    def ddt(self,t,x1x2):
        x1,x2 = np.hsplit(x1x2,2)
        x1 = x1.reshape(self.x1.shape)
        x2 = x2.reshape(self.x2.shape)
        A1x1 = self.A1.dot(x1)
        A2x2 = self.A2.dot(x2)
        if self.hybrid:
            x1tag = - x1 + self.beta[0] * self.apply_to_f(self.interaction(self.inter_term(x2,2),"interdependent")) * (1 - x1) *self.A1.dot(x1)
            x2tag = - x2 + self.beta[1] * self.apply_to_f(self.interaction(self.inter_term(x1,1),"competitive")) * (1 - x2) *self.A2.dot(x2)
            #x1tag = -x1 + self.beta * self.inter_term(x2,2) * (1 - x1) * self.A1.dot(x1)
            #x2tag = -x2 + self.beta * (1 - self.inter_term(x1,1)) * (1 - x2) * self.A2.dot(x2)

        else:
            if self.local_mf:
                x1tag = - x1 + self.beta[0] * self.apply_to_f(self.interaction(A2x2 / self.degvec2)) * (1 - x1) *A1x1
                x2tag = - x2 + self.beta[1] * self.apply_to_f(self.interaction(A1x1 / self.degvec1)) * (1 - x2) *A2x2
            else:
                x1tag = - x1 + self.beta[0] * self.apply_to_f(self.interaction(x2)) * (1 - x1) *A1x1
                x2tag = - x2 + self.beta[1] * self.apply_to_f(self.interaction(x1)) * (1 - x2) *A2x2
        if self.drivers:
            for i in range(self.N):
                if self.x1driver[i] == 1:
                    x1tag[i] = 0
                if self.x2driver[i] == 1:
                    x2tag[i] = 0
        return np.vstack((x1tag,x2tag)).flatten()

    def integrate(self,tfinal,steps=None):
        from scipy.integrate import solve_ivp
        if not steps:
            steps = tfinal / 0.01
        if self.diffusion_mode:
            ddt = self.diff_ddt
        else:
            ddt = self.ddt
        output = solve_ivp(ddt, (0, tfinal), np.vstack((self.x1, self.x2)).flatten(),method=self.solver)
        #x1x2,status = odeint(ddt,y0=np.vstack((self.x1,self.x2)).flatten(),t=np.linspace(0,tfinal,steps),full_output=True)
        self.status = output["status"]
        return output["y"]

    def integrate_by_chunks(self,tfinal,chunk_length=0.5,stop_condition=False):
        t=0
        op1 = np.array([])
        op2 = np.array([])
        while t < tfinal:
            x1x2 = self.integrate(chunk_length).T
            x1, x2 = np.hsplit(x1x2[1:], 2)
            self.x1 = x1[-1].reshape(self.x1.shape)
            self.x2 = x2[-1].reshape(self.x2.shape)
            op1 = np.concatenate((op1,self.order_param(x1.T,1)))
            op2 = np.concatenate((op2, self.order_param(x2.T, 2)))
            t+=chunk_length
            if stop_condition and (op1[-1] < 1e-4 or op2[-1] < 1e-3):
                break
        return np.vstack((op1,op2))

    def integrated_segments(self,nsegments=3,segment_time=2):
        segments=[]

        for segment in range(nsegments):
            x1x2 = self.integrate(tfinal=segment_time)
            x1, x2 = np.hsplit(x1x2, 2)
            self.x1 = x1[-1].reshape(self.x1.shape)
            self.x2 = x2[-1].reshape(self.x2.shape)
            segments.append([self.order_param(x1.T,1),self.order_param(x2.T,2)])
        return segments

    def betapath(self,betavec,tfinal):
        history=[]
        for beta in np.hstack((betavec,betavec[-2:0:-1])):
            self.beta = [beta,beta]
            x1x2 = self.integrate(tfinal,steps=None)
            x1, x2 = np.hsplit(x1x2[-1], 2)
            self.x1 = x1.reshape(self.x1.shape)
            self.x2 = x2.reshape(self.x2.shape)
            self.x1 += 0.001*np.random.random(self.x1.shape)
            self.x2 += 0.001*np.random.random(self.x2.shape)
            print("%.4f : (%.4f, %.4f)"%(beta,self.order_param(self.x1,1),self.order_param(self.x2,2)))
            history.append([self.order_param(self.x1,1),self.order_param(self.x2,2)])
        return history

    def betapath_from_points(self,p0,p1,steps,tfinal):
        p0 = np.array(p0); p1 = np.array(p1)
        m = p1 - p0
        history = []
        beta2vec = p0 + np.outer(np.linspace(0,1,steps),m)
        betavec2 = np.vstack((beta2vec,beta2vec[-2:0:-1]))
        oldop1, oldop2 = [self.order_param(self.x1, 1), self.order_param(self.x2, 2)]
        for beta2 in betavec2:
            self.beta = beta2
            #x1x2 = self.integrate(tfinal, steps=None)
            op1,op2 = self.integrate_by_chunks(tfinal)
            self.x1 += 0.001*np.random.random(self.x1.shape)
            self.x2 += 0.001*np.random.random(self.x2.shape)
            op1=op1[-1];op2=op2[-1]
            print("(%.6f, %.6f ) : (%.9f, %.9f) -> (%.9f, %.9f)" % (beta2[0],beta2[1], oldop1,oldop2,op1,op2))
            oldop1,oldop2 = [self.order_param(self.x1,1), self.order_param(self.x2,2)]
            history.append(np.array([op1,op2]).flatten())
        return betavec2,np.array(history)

    def betapath_from_points_global_op_compare(self,p0,p1,steps,tfinal):
        p0 = np.array(p0); p1 = np.array(p1)
        m = p1 - p0
        history = []
        global_op_history = []
        beta2vec = p0 + np.outer(np.linspace(0,1,steps),m)
        betavec2 = np.vstack((beta2vec,beta2vec[-2:0:-1]))
        oldop1, oldop2 = [self.order_param(self.x1, 1), self.order_param(self.x2, 2)]
        for beta2 in betavec2:
            self.beta = beta2
            #x1x2 = self.integrate(tfinal, steps=None)
            op1,op2 = self.integrate_by_chunks(tfinal)
            self.x1 += 0.001*np.random.random(self.x1.shape)
            self.x2 += 0.001*np.random.random(self.x2.shape)
            op1=op1[-1];op2=op2[-1]
            print("(%.6f, %.6f ) : (%.9f, %.9f) -> (%.9f, %.9f)" % (beta2[0],beta2[1], oldop1,oldop2,op1,op2))
            oldop1,oldop2 = [self.order_param(self.x1,1), self.order_param(self.x2,2)]
            history.append(np.array([op1,op2]).flatten())
            global_op_history.append(np.array([self.order_param(self.x1,1,True),self.order_param(self.x2,2,True)]).flatten())
        return betavec2,np.array(history),global_op_history

    def update_adjacency_matrix(self):
        self.A1 = nx.to_scipy_sparse_matrix(self.G1)
        self.A2 = nx.to_scipy_sparse_matrix(self.G2)
        self.degvec1 = np.array(self.A1.sum(axis=0)).reshape(self.x1.shape)
        self.degvec2 = np.array(self.A2.sum(axis=0)).reshape(self.x2.shape)
        self.degsum1 = sum(self.degvec1.flatten())
        self.degsum2 = sum(self.degvec2.flatten())
        self.L1 = nx.laplacian_matrix(self.G1)
        self.L2 = nx.laplacian_matrix(self.G2)

    def diff_ddt(self,x1x2,t):
        x1,x2 = np.hsplit(x1x2,2)
        x1 = x1.reshape(self.x1.shape)
        x2 = x2.reshape(self.x2.shape)
        x1tag = (self.ddriver1 - self.beta[0] * self.apply_to_f(self.interaction(self.inter_term(x2, 2))) * self.L1.dot(x1))*(self.x1*(1 - self.x1))
        x2tag = (self.ddriver2 - self.beta[1] * self.apply_to_f(self.interaction(self.inter_term(x1, 1))) * self.L2.dot(x2))*(self.x2*(1 - self.x2))
        return np.vstack((x1tag, x2tag)).flatten()

    def percolate(self, driver_fraction=0.1, nsamples = 100,tfinal=10,forward=True):
        G1 = self.G1.copy()
        G2 = self.G2.copy()
        self.add_drivers(driver_fraction=driver_fraction)
        self.randinit_upto([0.1,0.1])
        self.diffusion_mode = False
        G1edges = list(G1.edges())
        G2edges = list(G2.edges())
        np.random.shuffle(G1edges)
        np.random.shuffle(G2edges)
        if forward:
            self.G1 = nx.Graph()
            self.G2 = nx.Graph()
            self.G1.add_nodes_from(G1.nodes())
            self.G2.add_nodes_from(G2.nodes())

        maxedges = max(len(G1edges),len(G2edges))
        edgesatonce = maxedges // nsamples
        i=0
        res1,res2=[[],[]]
        while i<maxedges:
            if forward:
                self.G1.add_edges_from(G1edges[i:i + edgesatonce])
                self.G2.add_edges_from(G2edges[i:i + edgesatonce])
            else:
                self.G1.remove_edges_from(G1edges[i:i + edgesatonce])
                self.G2.remove_edges_from(G2edges[i:i + edgesatonce])

            self.update_adjacency_matrix()
            i+=edgesatonce
            op1,op2 = self.integrate_by_chunks(tfinal=tfinal)
            op1=op1[-1];op2=op2[-1]
            k1 = self.degsum1 / self.N
            k2 = self.degsum2 / self.N
            res1.append([k1,op1])
            res2.append([k2,op2])
            print()
            print(res1[-1]+[self.x1.mean(), len(max(nx.connected_components(self.G1), key=len))/self.N],end="\t")
            print(res2[-1]+[self.x2.mean(), len(max(nx.connected_components(self.G2), key=len))/self.N])
        self.G1 = G1.copy()
        self.G2 = G2.copy()
        self.update_adjacency_matrix()
        return res1,res2
