from magnumnp import *
import torch
from torch import sin, cos

# initialize mesh
n  = (50, 50, 1)
dx = (5e-9, 5e-9, 5e-9)
mesh = Mesh(n, dx)
h = 300e-3/constants.mu_0
t_f = 2e-9
dt = 1e-12

solver = "odeint_adjoint"
method = "implicit_adams"
loop = True

epoch_num = 100

# initialize material
material = {
           "Ms": 8e5,
           "A": 1.3e-11,
           "alpha": 0.08,
           "K": 1e5,
           "K_axis": [0,0,1]
           }

state = State(mesh, material)

state.m = state.Constant([0,0,0])
state.m[:, :, :, 2] = 1
state.m[:2, :, :, 0] = 0.1
state.m[2:, :, :, 0] = -0.1
state.m.normalize()
write_vti(state.m, "data/m0.vti", state)

m_t = state.Constant([0,0,0])
m_t[:, :, :, 0] = 1
m_t[:, :, :, 2] = 1
m_t.normalize()

demag    = DemagField(state)
exchange = ExchangeField(state)
phi = state.Tensor(state.Tensor([0.]), requires_grad = True)
theta = state.Tensor([1.5708], requires_grad = True)

external = ExternalField(state, torch.stack([sin(theta)*cos(phi),
                                             sin(theta)*sin(phi),
                                             cos(theta)], dim=-1))

aniso = AnisotropyField(state)

llg = LLGSolver(state, [demag, exchange, external, aniso])
llg.step(1e-10)
my_loss = torch.nn.L1Loss()
L = my_loss(state.m, m_t)
L.backward()

print("dL/dphi:", phi.grad)
print("dL/dtheta:", theta.grad)


#options = {"method": "dopri5",
#           "t0": 0.,
#           "t1": t_f,
#           "h": dt,
#           "rtol": 1e-3,
#           "atol": 1e-4,
#           "print_neval": False,
#           "neval_max": 1e6,
#           "safety": None}
#
## target magnetization orientation
#m_t = torch.zeros_like(m0)
#m_t[:, :, :, 0] = 1
#m_t[:, :, :, 2] = 1
#m_t = m_t / torch.norm(m_t, dim=3, keepdim=True)
#write_vtr(m_t, filename + '/target_state', mesh)
#
## =============================================================================
## def my_loss(m_res, m_tar):
##     losses = m_res - m_tar
##     return losses.pow(2).sum()
## =============================================================================
#
#my_loss = torch.nn.L1Loss()
#
#demag    = DemagField(mesh, material).to(cuda)
#exchange = ExchangeField(mesh, material).to(cuda)
#external = ExternalField_nn(mesh, material, h).to(cuda)
#aniso    = AnisotropyField(mesh, material).to(cuda)
#    
#model = Network([demag, exchange, external, aniso], mesh, material)
#model = model.to(cuda)
#optimizer = torch.optim.Adam(external.parameters(), lr=0.05) #0.05
#loss_list = []
#phi_list = []
#theta_list = []
#
#if os.path.exists(cache_dir+f"/m_relaxed_{n[0]}_{n[1]}_{n[2]}.pt"):
#    m_relaxed = torch.load(cache_dir+f"/m_relaxed_{n[0]}_{n[1]}_{n[2]}.pt")
#else:
#    m_relaxed = model.relax(m0, [demag, exchange, aniso])
#    torch.save(m_relaxed, cache_dir+f"/m_relaxed_{n[0]}_{n[1]}_{n[2]}.pt")
#write_vtr(m_relaxed, filename + "/m_relaxed", mesh)
#
#for epoch in range(epoch_num):
#    print('Epoch '+str(epoch))
#    optimizer.zero_grad()
#   # with torch.autograd.set_detect_anomaly(True):
#    result = model(m_relaxed, t_f, dt, solver, method=method, options=options, loop=loop)
#    
#    t = torch.arange(0,t_f,dt)
#    with open(datadir + '/nn_'+str(epoch)+'.dat', 'w') as f:
#        for i in range(result.shape[0]):
#            m = result[i,:,:,:]
#                #if i % 10 == 0:
#                 #   write_vtr(m, "data/nn_m_%04d" % (i/10))
#            f.write("%g %g %g %g\n" % ((t[i],) + tuple(torch.mean(m, axis=(0,1,2)))))
#   
#    write_vtr(result[-1], epochevoldir + '/final_state_'+str(epoch))#[-1]
#    
#    loss = my_loss(result[-1], m_t)#[-1] 
#
#    loss_list.append(loss.item())
#    phi_list.append(external.phi.item()*180/constants.pi)
#    theta_list.append(external.theta.item()*180/constants.pi)
#        
#    magn_plot_file(datadir + '/nn_'+str(epoch)+'.dat', picdir)
#    plot_param([loss_list, phi_list, theta_list], ["Loss", 'phi (deg)', 'theta (deg)'], filename)
#
#    start_time = tm()
#    loss.backward()
#    print('\n'+"--- Backward: %.3f seconds ---" % (tm() - start_time))
#    optimizer.step()
#    sys.exit()
#
#write_pvd(epochevoldir, '/epoch_evol', torch.arange(epoch_num))
#magn_plot_file(filename+'/nn_.dat', filename)
