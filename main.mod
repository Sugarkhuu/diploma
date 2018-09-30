var
lambda_p, c_p, h_p, q_h, d_p, w_p, l_p                              %patient
lambda_i, c_i, h_i, s_i, w_i, l_i, b_i   
lambda_e, c_e, q_k, s_e, y_e, r_k, u, k_e, x, l_ep, l_ei, b_ee      %entreprenuers 
pi, pi_wp, pi_wi, PIW                                               %
epsilon_d, epsilon_be, epsilon_bh, epsilon_qk, epsilon_l, epsilon_y, epsilon_z, epsilon_h, m_e, a_e, m_i, eps_K_b  
k, i, c, y, j_r, y1, c_g, tax, c_pr,                                % market clearing  
R_b,r_bh, r_be, r_t, r_d, K_b, B, BH, BE, D, b_h, b_e, d_b, j_b, J_b, v_b, spr_b, weight_BH, weight_BE;  % banking

varexo  e_z  e_h  e_mi  e_me  e_ae  e_l   e_qk  e_d  e_y  e_bh  e_be  e_r_t  e_eps_K_b;

parameters hab_p  beta_p  hab_i  beta_i   hab_e  delta  csi1  csi2  alpha  mu  beta_e   k_i   iota_w   piss  k_be  phi
           k_p  k_w     k_kb  k_bh  k_d gamma_p gamma_i gamma_e gamma_b  phi_R  phi_pi  phi_y delta_b iota_p
           rho_e_z  rho_e_h  rho_mi  rho_me  rho_e_ae rho_e_l  rho_e_qk  rho_e_y  rho_e_bh  rho_e_be  rho_e_d   rho_eps_K_b 
           m_i_ss   m_e_ss  epsilon_d_ss  epsilon_be_ss  epsilon_bh_ss  epsilon_y_ss   epsilon_l_ss  h eps_d  eps_bh eps_be pi_share
           r_t_ss  r_d_ss  r_bh_ss  r_be_ss  r_k_ss   J  eps_b vb_ss  eta pcg2pcp;

hab_p  =0;
beta_p =0.99455; 
hab_i  =0;
beta_i =0.975;  
hab_e  =0;
delta =0.03838; 
alpha =0.34; 
mu  = 0.8;
beta_e = beta_i;
phi = 1.5;
vb_ss  = 0.09*1.02;
eta    = 0.20*1.250;
iota_w =0.2757; 
iota_p =0.1605 ;           
k_p   =28.6502;
k_kb  =11.0683;
k_bh  =10.0867 ;
k_be  =9.3638;
k_d   =3.5030;
k_i   =10.1822;
k_w   =99.8983;
gamma_p =1;
gamma_i =1;
gamma_e =1;
gamma_b =1;
phi_R  =0.75;
phi_pi  =1.9816;
phi_y   =0.3459;           
rho_e_z  =0.95;
rho_e_h  =0.95;
rho_mi  =0.95;
rho_me  =0.95;
rho_e_ae =0.95;
rho_e_l  =0;
rho_e_qk  =0.95;
rho_e_y  =0;
rho_e_bh  =0;
rho_e_be  =0;
rho_e_d =0;
rho_eps_K_b  = 0.90;
piss  = 1;                                                         %steady state gross inflation rate
eps_d = - 1.1;                                                     % elast. of subst. of deposits 
eps_be = 5.85;                                                     % elast. of subst. of loans to H
eps_bh=  4.85;                                                     % elast. of subst. of loans to E
m_i_ss =0.75;                                                      % steady state loan-to-value ratio impatient households
m_e_ss =0.40;                                                      % steady state loan-to-value ratio Entrepreneurs
epsilon_d_ss     = (eps_d/(eps_d - 1));                            % steady state  mark down interest rate on deposit   =0.5935
epsilon_be_ss    =(eps_be/(eps_be - 1));                           % steady state mark up on loans to E    =1.4717
epsilon_bh_ss    =(eps_bh/(eps_bh - 1));                           % steady state mark up on loans to H    =1.5587
epsilon_y_ss = 8.9;
epsilon_l_ss = 5; 
r_d_ss  = (piss/beta_p) -1 ;
r_t_ss  = (piss/beta_p -1) *(eps_d - 1)/eps_d;                     % steady state gross nominal interest rate
r_bh_ss = epsilon_bh*r_t_ss ;                                      % steady state interest rate on loans to E
r_be_ss = epsilon_be*r_t_ss ;                                      % steady state interest rate on loans to H
r_k_ss  = -(1-delta)-m_e_ss*(1-delta)*piss/beta_e*
               (1/(1+r_be_ss)-beta_e/piss)+1/beta_e;               % steady state rental rate of capital
eps_b        = 0.5*(eps_bh + eps_be);
pi_share    = 0.55;                   % share of bank profits payed out to patient households OK!!
delta_b = r_t_ss/vb_ss * (eps_d - eps_b + v_b*eps_d*(eps_b-1))/((eps_b-1)*(eps_d-1)); 
csi1  = r_k_ss;
csi2  =0.1*r_k_ss;
h =1;                                                              % fixed supply housing
J = 0.203;
pcg2pcp = 0.32; 

model;

exp(lambda_p) = (1 -hab_p) * exp(epsilon_z) / (exp(c_p) - hab_p * exp(c_p(-1)));                                %%%/* dL / d c_p */
J*(exp(epsilon_h)) / exp(h_p)  - exp(lambda_p) * exp(q_h) + beta_p * exp(lambda_p(+1)) * exp(q_h(+1)) =0;       %%%/* dL / d h_p */
exp(lambda_p) = beta_p *  exp(lambda_p(+1)) * (1 + exp(r_d))/exp(pi(+1));   %%%/* dL / d d_p */
exp(c_p) + exp(q_h) * ( exp(h_p) - exp(h_p(-1))) + exp(d_p) = exp(w_p) * exp(l_p) + (1 + exp(r_d(-1)))* exp(d_p(-1))/exp(pi) + exp(j_r) / gamma_p + (1-pi_share)*exp(j_b(-1))/(gamma_p*exp(pi));         %%%/* BC*/
exp(lambda_i) = (1 - hab_i) * exp(epsilon_z) / ( exp( c_i) - hab_i * exp( c_i(-1))) ;                           %%%/* dL / d c_i */
J*( exp(epsilon_h)) / exp(h_i) - exp(lambda_i) * exp(q_h) + beta_i * exp(lambda_i(+1)) * exp(q_h(+1)) + exp(s_i) * exp( m_i) * exp(q_h(+1)) * exp(pi(+1)) = 0 ;     %%%  /* dL / d h_i*/
exp(lambda_i)  -  beta_i *  exp( lambda_i ( +1)) *( 1 + exp ( r_bh)) / exp(pi(+1)) = exp(s_i) * ( 1 + exp(r_bh)) ;                           %%/* dL / d b_i */
exp(c_i) + exp(q_h) * (exp(h_i) - exp(h_i(-1))) + ( 1 + exp(r_bh(-1))) * exp(b_i(-1)) / exp(pi) = exp(w_i) *  exp(l_i) + exp(b_i) ;          %%%                                                                                
(1 + exp(r_bh)) * exp(b_i) = exp(m_i) * exp(q_h(+1)) * exp(h_i) * exp( pi(+1));      %%%/*BORROWING CONSTRAINT for IH*/
exp(lambda_e) = (1 - hab_e) / ( exp(c_e) - hab_e * exp(c_e(-1))) ;                                              %%/* dL / d c_e */
exp(lambda_e) * exp(q_k) = exp(s_e) * exp(m_e) * exp (q_k(+1)) * exp(pi(+1)) * (1 - delta) + beta_e * exp(lambda_e(+1)) * ( exp(r_k(+1)) * exp(u(+1)) + exp(q_k(+1))* (1 - delta) 
       - (csi1 * (exp(u(+1)) - 1) + csi2/2 *((exp(u(+1)) - 1)^2))) ;    %%%   /*dL / d k_e*/  
exp(r_k) = csi1 + csi2 *( exp(u) -1 ) ;                                                                         %%/* dL / d u   DEGREE OF UTILIZATION OF CAPITAL*/
exp(r_k) = alpha * exp(a_e) * exp(k_e(-1))^(alpha-1) * exp(u)^(alpha -1) * ( exp(l_ep)^mu * exp(l_ei)^(1-mu) ) ^ (1-alpha) / exp(x) ;  %%%/* MARGINAL PRODUCT OF CAPITAL */
exp(y_e) = exp(a_e) * ( exp(k_e(-1)) * exp(u))^(alpha) * ( exp(l_ep)^mu * exp(l_ei)^(1-mu) ) ^ (1-alpha) ;      %%%/* PRODUCTION FUNCTION*/
exp(w_p) = (1 - alpha) * exp(y_e) * mu / (exp(x) * exp(l_ep) );               %%%/* dL / d l_ep */
exp(w_i) = (1 - alpha) * exp(y_e) * (1 - mu) / (exp(x) * exp(l_ei));          %%%/* dL / d l_ei */
exp(lambda_e) - exp(s_e) *(1 + exp(r_be))  =  beta_e * exp(lambda_e(+1)) * ( 1 + exp(r_be)) / exp(pi(+1)) ;      %%%/* dL / d b_e */
exp(c_e) + (exp(w_p)*exp(l_ep) + exp(w_i)*exp(l_ei)) + ((1 + exp(r_be(-1))) * exp(b_ee(-1)) / exp(pi))  + exp(q_k) * exp(k_e) +
(csi1*(exp(u)  - 1) + csi2/2*(exp(u) -1)^2) * exp(k_e(-1)) = exp(y_e)/exp(x) + exp(b_ee) + exp(q_k) *(1 - delta) * exp(k_e(-1)) ;             %%%      /*BC*/
(1 + exp(r_be)) * exp(b_ee) = exp(m_e) * exp( q_k(+1)) * exp( pi(+1)) * (1 - delta) * exp(k_e)  ;       %%%/* BORROWING CONSTRAINT*/
(1 -exp(epsilon_l)) * exp(l_p)+ exp(l_p)^(1 + phi) / exp(w_p) *exp(epsilon_l)/ exp(lambda_p) -
k_w * (exp(pi_wp) - exp(pi(-1))^iota_w * piss^(1 - iota_w)) * exp(pi_wp) + beta_p * exp(lambda_p(+1)) / exp(lambda_p) * k_w * ( exp(pi_wp(+1))  - exp(pi)^( iota_w) * piss^( 1 - iota_w)) * exp(pi_wp(+1))^(2) / exp(pi) = 0;
exp(pi_wp) = exp(w_p) / exp(w_p(-1)) * exp(pi) ;                    %%% /* WAGE INFLATION FOR PATIENT HOUSEHOLDS*/
(1 -exp(epsilon_l)) * exp(l_i) + exp(l_i)^(1 + phi) / exp(w_i) *exp(epsilon_l) / exp(lambda_i) - 
k_w * (exp(pi_wi) - exp(pi(-1))^iota_w * piss^(1 - iota_w)) * exp(pi_wi) + beta_i *  exp(lambda_i(+1)) / exp(lambda_i) * k_w * ( exp(pi_wi(+1))  - exp(pi)^(iota_w) * piss^( 1 - iota_w)) * exp(pi_wi(+1))^(2) / exp(pi) = 0 ;
exp(pi_wi) = exp(w_i) / exp(w_i(-1)) * exp(pi) ;            %%/* WAGE INFLATION FOR IMPATIENT HOUSEHOLDS*/
exp(k) = (1 - delta) * exp(k(-1)) + ( 1 - k_i/2 *(exp(i) * exp(epsilon_qk)/exp(i(-1)) - 1 )^2)* exp(i);            %%%/*AMOUNT OF NEW CAPITAL THAT CGP FIRMS CAN PRODUCE */
exp(q_k) * (1 - k_i/2 * (exp(i) * exp(epsilon_qk)/ exp(i(-1)) - 1) ^2 - 
k_i*(exp(i) * exp(epsilon_qk)/exp(i(-1)) - 1) * exp(i) *exp(epsilon_qk)/exp(i(-1))) + 
beta_e *  exp(lambda_e(+1))/exp(lambda_e) * exp(q_k(+1)) * k_i * ( exp(i(+1)) * exp(epsilon_qk(+1)) / exp(i)  - 1 )* exp(epsilon_qk(+1)) * ( exp(i(+1)) / exp(i))^2 = 1;  %%% /*REAL PRICE OF CAPITAL*/
1 - exp(epsilon_y) + 
exp(epsilon_y) / exp(x) - 
k_p * (exp(pi) -      (exp(pi(-1))^iota_p * piss^(1 - iota_p) ))* exp(pi) +
beta_p * (exp(lambda_p(+1))/ exp(lambda_p)) *  k_p * (exp(pi(+1)) - (exp(pi)^     iota_p* piss^(1 - iota_p)))* exp(pi(+1)) *(exp(y(+1)) / exp(y)) = 0 ;    %%% /* NON LINEAR PHILLIPS CURVE*/
exp(j_r) = exp(y) *(1 - (1/ exp(x)) - (k_p/2)*(exp(pi) - (exp(pi(-1))^(iota_p) * piss^(1 - iota_p)))^2);           %%%/* RETAILERS PROFITS*/
exp(R_b) = exp(r_t) - k_kb * ( (exp(K_b) / (exp(B)))  - exp(v_b) ) *( exp(K_b) / (exp(B)))^2;     %%%/*FOCS FROM THE MAXIMIZATION PROBLEM OF THE WHOLESALE BRANCH*/
/*exp(r_t) = exp(R_d) ;  *//*NO ARBITRAGE CONDITION WITH THE LENDING FACILITY OF THE CENTRAL BANK*/

exp(v_b) =  (vb_ss) ^( 1-0.92 )  *  ( exp(y1) / exp(y1(-4) ))^ ( 1 - 0.92 )  *  exp( v_b(-1) )^0.92;
exp(weight_BH) = ( 1 ^( 1-0.94 ) ) * ( exp(y1) / exp(y1(-4)) )^( ( 1 - 0.94 ) * ( - 10 * 1 ) ) * exp( weight_BH(-1) )^0.94; //32
exp(weight_BE) = ( 1 ^( 1-0.92 ) ) * ( exp(y1)/exp(y1(-4)) )^(( 1 - 0.92 ) * ( - 15 * 1  ))* exp( weight_BE(-1) )^0.92; //33

1 - exp(epsilon_bh)/(exp(epsilon_bh) -1)  + exp(epsilon_bh)/(exp(epsilon_bh) -1)  * exp(R_b) / exp(r_bh) - k_bh*(exp(r_bh)/exp(r_bh(-1))  -1 ) * exp(r_bh)/exp(r_bh(-1))  +  
beta_p * (exp(lambda_p(+1)) / exp(lambda_p)) * k_bh * ( exp(r_bh(+1)) / exp(r_bh)  - 1 )  *(( exp(r_bh(+1)) / exp ( r_bh))^2) * (exp(b_h(+1))/ exp(b_h)) = 0 ; 

1 - exp(epsilon_be)/(exp(epsilon_be) -1)   +  exp(epsilon_be)/(exp(epsilon_be) -1)  * exp(R_b) / exp(r_be) - k_be*(exp(r_be)/exp(r_be(-1))  -1 ) * exp(r_be)/exp(r_be(-1))  +  
beta_p * (exp(lambda_p(+1)) / exp(lambda_p)) * k_be * ( exp(r_be(+1)) / exp(r_be)  - 1 )  *(( exp(r_be(+1)) / exp ( r_be))^2) * (exp(b_e(+1))/ exp(b_e)) = 0 ; 

- 1 + exp(epsilon_d)/(exp(epsilon_d) -1)   -  exp(epsilon_d)/(exp(epsilon_d) -1)  * exp(r_t) / exp(r_d) - k_d*(exp(r_d)/exp(r_d(-1))  -1 )  * exp(r_d)/exp(r_d(-1))  +  
beta_p * (exp(lambda_p(+1)) / exp(lambda_p)) * k_d * ( (exp(r_d(+1)) / exp(r_d))  - 1 )   *(( exp(r_d(+1)) / exp ( r_d))^2)   * (exp(d_b(+1))/ exp(d_b)) = 0 ;   

exp(j_b)  = + exp(r_bh) * exp(b_h) 
            + exp(r_be) * exp(b_e)
            + exp(r_t)  * exp(B)*eta
            - exp(r_d)  * exp(d_b) 
            - exp(K_b)- k_kb/2 * ( ((exp(K_b) / (exp(B)))  - exp(v_b) ) ^2) * exp(K_b);                    
exp(K_b) * exp(pi) = (1 - delta_b) * exp(K_b(-1))  + (pi_share)*exp(j_b(-1)) ;             %%%          /*BANK CAPITAL ACCUMULATION*/
gamma_b * exp(d_b)  = gamma_p * exp(d_p) ;                %%%                      
gamma_b * exp(b_h)  = gamma_i * exp(b_i) ;                %%%
gamma_b * exp(b_e)  = gamma_e * exp(b_ee);                %%%
exp(B)*(1 - eta)  =  exp(d_b) + exp(K_b)  + eps_K_b;       %%%/*BALANCE SHEET CONSTRAINT  B = D + K_b*/
(1 + exp(r_t)) = (1 + r_t_ss)^(1 - phi_R) * (1 + exp(r_t(-1)))^(phi_R) *
                 ((exp(pi)/piss)^(phi_pi) * (exp(y1)/exp(y1(-1)))^phi_y)^(1 - phi_R) *(1+e_r_t) ;    %%%

//exp(y)  = exp(c) + 
//          /*exp(q_k)* */    ( exp(k) - (1 - delta) *exp(k(-1))) ;
//          /*+ exp(k_e(-1)) * (csi1*(exp(u) -1 ) + (csi2/2)*(exp(u) -1)^(2)) + */
//          /*delta_b* exp(K_b(-1)) / exp(pi)  */
         

//          /*- (k_p/2) * (exp(pi) /   - exp(pi(-1))^iota_p * piss^(1-iota_p)) ^(2)*exp(y) */
//          /*- k_d/2  * ( (exp(r_d)/exp(r_d(-1))-1)^2)   * exp(r_d) *exp(d_b)  */
//          /*- k_be/2 * ( (exp(r_be)/exp(r_be(-1))-1)^2) * exp(r_be)*exp(b_e)  */
//          /*- k_bh/2 * ( (exp(r_bh)/exp(r_bh(-1))-1)^2) * exp(r_bh)*exp(b_h);  */ 
exp(c_pr)= gamma_p * exp(c_p) + gamma_i * exp(c_i) + gamma_e * exp(c_e);
exp(c_g) = pcg2pcp*exp(c_p);
exp(c)   = gamma_p * exp(c_p) + gamma_i * exp(c_i) + gamma_e * exp(c_e) + exp(c_g)  ;%%
h        = gamma_p * exp(h_p) + gamma_i * exp(h_i) ;           %%
exp(k)   = gamma_e * exp(k_e) ;
exp(B)*(1 - eta)   = (exp(BH) + exp(BE)) ;
exp(BH)  = gamma_b * exp(b_h);
exp(BE)  = gamma_b * exp(b_e);
exp(D)   = gamma_p * exp(d_p);
exp(y)   = gamma_e * exp(y_e) ;
exp(J_b) = gamma_b * exp(j_b);  
gamma_e * exp(l_ep) = gamma_p * exp(l_p);
gamma_e * exp(l_ei) = gamma_i * exp(l_i);

exp(tax) = exp(c_g) + (eta)*( exp(d_b(-1)-exp(d_b) - exp(K_b)) + exp(K_b(-1)) )*(exp(r_t(-1))) ; //gov.equation
exp(spr_b)    = 0.5*exp(r_bh) + 0.5*exp(r_be) - exp(r_d); //Int.rate spread

exp(PIW)            = ( exp(w_p) + exp(w_i) )  / ( exp(w_p(-1)) + exp(w_i(-1)) ) * exp(pi);
exp(y1)             = exp(c) +    1     * (exp(k)-(1-delta)*exp(k(-1))) ;
                   //13

/***************************************************************+ EXOGENOUS PROCESS*****************************************************************************/

exp(epsilon_z)     = 1 - rho_e_z   *    1               + rho_e_z   * exp(epsilon_z(-1))    + e_z;        %%
exp(a_e)           = 1 - rho_e_ae  *    1               + rho_e_ae  * exp(a_e(-1))          + e_ae;       %%
exp(epsilon_h)     = 1 - rho_e_h   *    1               + rho_e_h   * exp(epsilon_h(-1))    + e_h;        %%
exp(m_i)           = (1-rho_mi)    *  m_i_ss            + rho_mi    * exp(m_i(-1))          + e_mi;       %%
exp(m_e)           = (1-rho_me)    *  m_e_ss            + rho_me    * exp(m_e(-1))          + e_me;       %%
exp(epsilon_d)     = (1-rho_e_d)   * epsilon_d_ss       + rho_e_d   * exp(epsilon_d(-1))    + e_d;        %%
exp(epsilon_bh)    = (1-rho_e_be)  * epsilon_be_ss      + rho_e_be  * exp(epsilon_be(-1))   + e_be;       %%
exp(epsilon_be)    = (1-rho_e_bh)  * epsilon_bh_ss      + rho_e_bh  * exp(epsilon_bh(-1))   + e_bh;       %%
exp(epsilon_qk)    =  1-rho_e_qk   *    1               + rho_e_qk  * exp(epsilon_qk(-1))   + e_qk;       %%
exp(epsilon_y)     = (1-rho_e_y)   * epsilon_y_ss       + rho_e_y   * exp(epsilon_y(-1))    + e_y;        %%
exp(epsilon_l)     = (1-rho_e_l)   * epsilon_l_ss       + rho_e_l   * exp(epsilon_l(-1))    + e_l;        %%
exp(eps_K_b)       = (1-rho_eps_K_b)*    1              + rho_eps_K_b* exp(eps_K_b(-1))     + e_eps_K_b;  %%
end;

initval;
lambda_p=-0.17442470738849;
c_p=0.174368372837109;
h_p=-0.0615335201724342;
q_h=3.85353808148559;
d_p=1.48638779662852;
w_p=0.121577429545804;
l_p=-0.184040916717857;
lambda_i=1.56059206530034;
c_i=-1.5605733973416;
h_i=-2.81821787495264;
s_i=-2.76908841457159;
w_i=-1.40415122176564;
l_i=-0.0446113806026163;
b_i=0.735920446387513;
lambda_e=2.31067950236279;
c_e=-2.31068210318773;
q_k=8.97518352660628E-07;
s_e=-2.05353647417474;
y_e=0.695380284527115;
r_k=-2.85118627163358;
u=0.516679968787456;
k_e=1.83166779489747;
x=0.119196587211381;
l_ep=-0.184040128082612;
l_ei=-0.044615149986417;
b_ee=0.863993120842196;
pi=-0.000477094131243953;
pi_wp=-0.000479950332114822;
pi_wi=-0.000480165469969359;
PIW=-0.000480007225928979;
epsilon_d=-0.646627107399912;
epsilon_be=0.230905604112095;
epsilon_bh=0.187463923179717;
epsilon_qk=6.20206955585757E-09;
epsilon_l=1.60943791248643;
epsilon_y=2.18605127685189;
epsilon_z=-5.34448976957994E-09;
epsilon_h=0.000011647006010751;
m_e=-0.916287639617494;
a_e=1.07225774443232E-06;
m_i=-0.287676963434576;
eps_K_b=-3.17270665061223E-08;
k=1.83168851626454;
i=-1.42884759284115;
c=0.631614407266322;
y=0.695379729284502;
j_r=-1.49016722472392;
y1=0.751095838599898;
c_g=-0.964823909226989;
tax=-0.964644315979603;
c_pr=0.405285562853627;
R_b=-4.63810571050966;
r_bh=-4.45064190699467;
r_be=-4.40720004382034;
r_t=-4.64397751396188;
r_d=-5.29060463107096;
K_b=-3.23518557998594;
B=1.78249662124258;
BH=0.735995645207136;
BE=0.863852713680075;
D=1.48640827406987;
b_h=0.735951621668061;
b_e=0.863963300857266;
d_b=1.48635003753294;
j_b=-5.1901304661428;
J_b=-5.19010870211428;
v_b=-2.38814298128211;
spr_b=-4.97622537717927;
weight_BH=-5.1520256177705E-17;
weight_BE=6.50777245928301E-16;

end;

resid(1);
//steady(solve_algo=0);
check;




          % e_z   % e_h   % e_mi  % e_me  % e_a  % e_l    % e_qk   % e_y% %% e_bh     % e_be   %e_d    %e_r_t    %e_eps_K_b 
Std_e = [0.027     0      0       0       0       0        0        0      0         0        0          0          0              % e_z
          0      0.076    0       0       0       0        0        0      0         0        0          0          0              % e_h
          0        0     0.003    0       0       0        0        0      0         0        0          0          0              % e_mi
          0        0      0      0.007    0       0        0        0      0         0        0          0          0              % e_me
          0        0      0       0      0.06     0        0        0      0         0        0          0          0              % e_a
          0        0      0       0       0     0.57       0        0      0         0        0          0          0              % e_l
          0        0      0       0       0       0       0.019     0      0         0        0          0          0              % e_qk
          0        0      0       0       0       0        0      0.634    0         0        0          0          0              % e_y
          0        0      0       0       0       0        0        0     0.067      0        0          0          0              % e_bh 
          0        0      0       0       0       0        0        0      0         0.063    0          0          0              % e_be
          0        0      0       0       0       0        0        0      0         0        0.033      0          0              %e_d
          0        0      0       0       0       0        0        0      0         0        0          0.002      0              %e_r_t
          0        0      0       0       0       0        0        0      0         0        0          0          0.031];        %e_eps_K_b
                


shocks;
var e_z   =  Std_e(1,1)^(2);
var e_h   =  Std_e(2,2)^(2);
var e_mi  =  Std_e(3,3)^(2);
var e_me  =  Std_e(4,4)^(2);
var e_ae  =  Std_e(5,5)^(2);
var e_l   =  Std_e(6,6)^(2);
var e_qk  =  Std_e(7,7)^(2);
var e_y   =  Std_e(8,8)^(2);
var e_bh  =  Std_e(9,9)^(2);
var e_be  =  Std_e(10,10)^(2);
var e_d   =  Std_e(11,11)^(2);
var e_r_t =  Std_e(12,12)^(2);
var e_eps_K_b =Std_e(13,13)^(2);
end;


options_.nograph   = 0;
options_.nomoments = 0;
options_.noprint   = 0;

stoch_simul(order=1,irf=40) y1 y c i B q_h;

%% stoch_simul(order=1,irf=40) y1 y c i B q_h q_k J_b BE BH K_b  pi  r_bh r_be r_d r_t D  h_i h_p v_b;


/*
ss_vector = oo_.steady_state;
n_vars    = length(ss_vector);
ss_cell   = mat2cell(ss_vector,ones(1,n_vars));
tmp       = [M_.endo_names,ones(n_vars,1)*44]';
varnames  = ['[',tmp(:)',']'];
eval([varnames,' = deal(ss_cell{:});'])


disp(' ');disp('%%% Display some SS results : %%%%')

disp(['C/Y:            ',num2str(exp(c)/exp(y1))]);
disp(['C_Private/Y:    ',num2str(exp(c_pr)/exp(y1))]);
disp(['G/Y:            ',num2str(exp(c_g)/exp(y1))]);
disp(['I/Y:            ',num2str(exp(i)/exp(y1))]);
disp(['K/Y:            ',num2str(exp(k)/exp(y1))]);
disp(['R_b  (%annual): ',num2str(exp(R_b)*400)]);
disp(['r_t  (%annual): ',num2str(exp(r_t)*400)]);
disp(['r_bh (%annual): ',num2str(exp(r_bh)*400)]);
disp(['r_be (%annual): ',num2str(exp(r_be)*400)]);
disp(['r_d  (%annual): ',num2str(exp(r_d)*400)]);
disp(['spread(%annual):',num2str(exp(spr_b)*400)]);
disp(['B:              ',num2str(exp(B))]);
disp(['B/Y:            ',num2str(exp(B)/exp(y1))]);
disp(['D/Y:            ',num2str(exp(d_p)/exp(y1))]);
disp(['K_b/B:          ',num2str(exp(K_b)/exp(B))]);
disp(['K_b/Y:          ',num2str(exp(K_b)/exp(y1))]);
disp(['Y:              ',num2str(exp(y1))]);

steady_state_variables(1,:) = [ exp(y1),... ; % 1
    exp(c),...                           ; % 2
    exp(i),...                           ; % 3
    exp(c_p),...                         ; % 4
    exp(c_i),...                         ; % 5
    exp(c_e),...                         ; % 6
    exp(h_p),...                         ; % 7
    exp(h_i),...                         ; % 8
    exp(l_p),...                         ; % 9
    exp(l_i),...                         ; % 10
    exp(BE),...                          ; % 11
    exp(BH),...                          ; % 12
    exp(q_h),...                         ; % 13
    exp(K_b),...                         ; % 14
    exp(y),...                           ; % 15
    exp(spr_b),...                       ; % 16
    exp(B),...                           ; % 17
    exp(D),...                           ; % 18
    M_.params(30),...                    ; % 19
    exp(J_b)]                            ; % 20

volatilities(1,:) = sqrt(diag(oo_.var))'.*100;

for ii = 1:length(oo_.var);
    for jj = 1:length(oo_.var);
        corr_mat(ii,jj) = oo_.var(ii,jj)/( sqrt(oo_.var(ii,ii)) * sqrt(oo_.var(jj,jj)) );
    end
end

welfare =       ( 1/( 1-beta_p ) )* ( log( exp(c_p)*( 1-hab_p ) ) + j * log( exp(h_p) ) - ( 1 / ( 1 + phi ) ) * exp(l_p)^( 1 + phi ) ) +...
                ( 1/( 1-beta_i ) )* ( log( exp(c_i)*( 1-hab_p ) ) + j * log( exp(h_i) ) - ( 1 / ( 1 + phi ) ) * exp(l_i)^( 1 + phi ) ) +...
                ( 1/( 1-beta_e ) )* ( log( exp(c_e)*( 1-hab_p ) ) );

*/
