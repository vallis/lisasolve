from mpi4py import MPI
from countdown import countdown

class masterslave(object):
    def __init__(self,inputlen,outputlen,countdown=True):
        self.size, self.rank = MPI.COMM_WORLD.Get_size(), MPI.COMM_WORLD.Get_rank()
        self.slaves = size - 1
        
        if self.rank == 0:
            self.send_bufs = [N.zeros(inputlen,'d') for s in range(self.slaves)]
            self.send_reqs = [MPI.COMM_WORLD.Send_init((send_bufs[s], MPI.DOUBLE), s+1, 0) for s in range(self.slaves)]
            
            self.recv_bufs = [N.zeros(outputlen,'d') for s in range(self.slaves)]
            self.recv_reqs = [MPI.COMM_WORLD.Recv_init((recv_bufs[s], MPI.DOUBLE), s+1, 1) for s in range(self.slaves)]
            
            self.countdown = countdown
            self.busy = [-1 for s in range(self.slaves)]
        else:
            self.recv_buf = N.zeros(inputlen,'d')            
            self.recv_req = MPI.COMM_WORLD.Recv_init((recv_buf, MPI.DOUBLE), 0, 0)
            
            self.send_buf = N.zeros(outputlen,'d')
            self.send_req = MPI.COMM_WORLD.Send_init((send_buf, MPI.DOUBLE), 0, 1)
    
    def itermaster(self,data,function):
        if self.rank != 0: return
        
        cnt = 0; ln = len(data)
        if self.countdown: c = countdown(ln,10000)
        while cnt < ln:
            for s in range(self.slaves):
                if self.busy[s] == -1 and cnt < ln:
                    self.send_bufs[s][:] = data[cnt][:]
                    self.send_reqs[s].Start()
                    self.recv_reqs[s].Start()
                    
                    self.busy[s] = cnt
                    cnt = cnt + 1
            
            if sum(x == -1 for x in self.busy) != len(self.busy):     # at least one of busy[s] is not -1
                f = MPI.Prequest.Waitany(self.recv_reqs)
                
                function(cnt,self.recv_bufs[f])
                
                busy[f] = -1
                c.status()
        c.end()
    
    def iterslave(self,function):
        if self.rank == 0: return
        
        while True:
            self.recv_req.Start()
            MPI.Prequest.Wait(self.recv_req)
            
            if recv_buf[0] == 0.0:
                return
            else:
                send_buf[:] = function(recv_buf)
                send_req.Start()
    
    def __del__(self):
        if self.rank == 0:
            for r in self.recv_reqs: r.Free()
            for r in self.send_reqs: r.Free()
        else:
            self.recv_req.Free(); self.send_req.Free()        

