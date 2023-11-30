# Cosmological Dynamic Systems
A collection of cosmological dynamical systems, simulated in MatLab. These were created to investigate a discrepancy in the number of critical points found in two different dynamical systems that were derived from the same fundamental cosmological equations.

## Preliminaries
In 1915 Einstein introduced his field equations 

$$\begin{equation}
    G_{\mu \nu} \equiv R_{\mu \nu} - \frac{1}{2}g_{\mu \nu}R =8\pi GT_{\mu \nu}
\end{equation}$$

where $G_{\mu \nu}$ is the Einstein tensor, $R_{\mu \nu}$ is the Ricci tensor, and $R$ is the Ricci scalar.  If we work in natural units ($G=c=1)$ and add the cosmological constant $\Lambda$ we get our modified field equation 

$$\begin{equation}R_{\mu \nu}-\frac{1}{2}g_{\mu \nu}R+\Lambda g_{\mu \nu} = 8\pi T_{\mu \nu}\end{equation}$$

Measurements of the cosmic microwave background radiation by NASA's Wilikinson Microwave Anisotropy Probe (WAMP) have reported the temperature at different parts of the universe all seem to be approximately 2.7260 kelvin, this means that on large enough scales (>100 MPc) the universe appears look the same in all directions.  

This data justifies us in assuming that the universe is homogeneous and isotropic on large scales, this is the *cosmological principle*.  Homogeneity implies that the universe has the same expansion rate everywhere while being isotropic implies the universe looks the same from all directions i.e. is invariant under rotation.  We also assume the following axioms:

1. The laws of physics do not change and are the same everywhere
2. Physical constants are true constants
3. The universe is connected

In this thesis we will consider a universe that assumes the cosmological principal which is called the Friedmann-Lemaître-Robertson-Walker (FLRW) universe with the line element $$ds^2 = -dt^2+a^2(t)[d\chi^2 +\Sigma(\chi)^2 d\Omega^2 ]$$

We are particularly interested in the case where the universe to be spatially flat ($k=0$).  We can then define the *Hubble parameter* 

$$H = \frac{\dot{a}}{a}$$
and then combining (1) and (2) we obtain the Friedmann equations $$\begin{align}
    H^2 &= \frac{8\pi\rho}{3} +\frac{\Lambda}{3}\\
    \frac{\ddot{a}}{a} &= -\frac{4\pi}{3}(\rho+3P)+\frac{\Lambda}{3}
\end{align}$$

while the continunity equation is given by $$\begin{equation}
    \dot{\rho}+ 3H(\rho + P) = 0
\end{equation}$$

where $\rho$ and $P$ are the total energy density and pressure of the fluid respectively, we can relate these two values though our equation of state $$\begin{equation}
    P = \rho w
\end{equation}$$ When our radiation and matter densities are uncoupled our continuity equation can be split in two $$\begin{align}
\dot{\rho_m} + 3H(\rho_m + P_m) &= 0\\
\dot{\rho_{rad}} + 3H(\rho_{rad} + P_{rad}) &= 0
\end{align}$$

With these fundemental equations we can begin to derive our first cosmological dynamical system.

## Cosmological dynamical system


Beginning with one of our Friedman equations with matter ($P_m = 0$) and radiation $(P_{rad} = \frac{1}{3}\rho_{rad})$:

$$3H^2 = \kappa^2\rho_m + \kappa^2\rho_{rad} + \Lambda$$ 
then we divide through by $3H^2$ to get our dimensionless variables $$\begin{equation}1 = \frac{\kappa^2\rho_m}{3H^2} +\frac{\kappa^2\rho_{rad}}{3H^2} + \frac{\Lambda}{3H^2}\end{equation}$$ $$\begin{align}
    x &= \frac{\kappa^2\rho_m}{3H^2} \\
    y &= \frac{\kappa^2\rho_{rad}}{3H^2} \\
    \Omega_{\Lambda} &= \frac{\Lambda}{3H^2} = \frac{\kappa^2\rho_{\Lambda}}{3H^2}
\end{align}$$where $\rho_{\Lambda} = \Lambda / \kappa^2$.  We can see that $x,y \geq 0$.  Now our field equation has become $$1 = x+ y + \Omega_{\Lambda} $$ Note that we can always substitute $\Omega_{\Lambda}$ for $x$ or $y$ so hence we are working in 2D system.  From this equation we also get the constraint $x+y \leq 1$ and so our state space forms a triangle $$\textbf{S} = \{(x,y) | 0\leq x,y  \text{ and } x+y \leq 1\}$$

To construct our first cosmological dynamical system we differentiate our variables with respect to time $t$ 
$$\begin{align*}
    \dot{x} &= \frac{\kappa^2\dot{\rho_m}}{3H^2} - \frac{2\kappa^2\rho_m}{3H}\frac{\dot{H}}{H}\\
    \dot{y} &= \frac{\kappa^2\dot{\rho_{rad}}}{3H^2}-\frac{2\kappa^2\rho_{rad}}{3H}\frac{\dot{H}}{H}
\end{align*}$$
Then we transform into differentiating by conformal time $\eta$ e to get our final dynamical system 

$$\begin{align*}
    x' = \frac{\text{d}x}{\text{d}\eta} = \frac{1}{H}\frac{\text{d}x}{\text{d}t} &= \frac{\kappa^2\dot{\rho_m}}{3H^3} - \frac{2\kappa^2\rho_m}{3H^2}\frac{\dot{H}}{H}\\
    y' = \frac{\text{y}x}{\text{d}\eta} = \frac{1}{H}\frac{\text{d}y}{\text{d}t} &= \frac{\kappa^2\dot{\rho_{rad}}}{3H^3}-\frac{2\kappa^2\rho_{rad}}{3H^2}\frac{\dot{H}}{H}
\end{align*}$$

To develop these equations we can use our conservation equation.  If we assume that matter and radiation are non-interacting then we get the following 
$$\begin{align*}
    \dot{\rho_m} + 3H(\rho_m + P_m) = \dot{\rho_m} + 3H\rho_m\phantom{...} &\rightarrow\phantom{...} \dot{\rho_m} = -3H\rho_m \\
    \dot{\rho_{rad}} + 3H(\rho_{rad} + \frac{1}{3}\rho_{rad}) = \dot{\rho_{rad}} + 3H(\frac{4}{3}\rho_{rad})\phantom{...} &\rightarrow \phantom{...}\dot{\rho_{rad}} = -4H\rho_{rad}
\end{align*}$$
We can use our acceleration equation to get an expression for $\frac{\dot{H}}{H}$ 
$$\frac{\dot{H}}{H} = \frac{-\kappa p^2}{2H^2} - \frac{3}{2} - \frac{\kappa^2}{2a^2H^2} = -\frac{1}{2}H(3-y-\Omega_{\Lambda})$$

Subbing this all in gives our final 2D dynamical system 
$$\begin{align*}
    x' &= x(3x+4y-3)\\
    y' &= y(3x+4y-4)
\end{align*}$$
This system has three critical points: $$\begin{align}
R &= (0,1), \text{ a radiation dominated universe}\\
O &= (0,0), \text{ a dark-energy dominated universe}\\
M &= (1,0), \text{ a matter dominated universe}
\end{align}$$
with Jacobian $$\begin{pmatrix}
6x+4y-3 & 4x \\
3y & 3x+8y-4
\end{pmatrix} $$Now we can perform linear stability analysis and evaluate the Jacobian at each critical point to understand their respective stability. $$J_R = \begin{pmatrix}
3& 4 \\
0 & -1
\end{pmatrix} \phantom{.....}J_O = \begin{pmatrix}
-3& 0 \\
0 & -4
\end{pmatrix}\phantom{.....}J_M = \begin{pmatrix}
1& 0 \\
3 &8
\end{pmatrix}$$These give the following eigenvalues: $R: \{-1,3\}$, $O: \{-3,-4\}$, and $M: \{1,4\}$, hence we can classify $R$ as saddle point, $O$ as a stable node, and $M$ as an unstable node.  We can see the phase portrait generated by Matlab of this system

![alt text](https://github.com/BillGriffin98/CosmologicalDynamicSystems/blob/main/orginaltrianglepicture.png)

From this portrait we can see the hetroclinic trajectory $R \rightarrow M \rightarrow O$ which describes a universe that begins in a radiation dominated era, progresses to a matter dominated era, before finally tending towards a dark-matter dominated universe.  A trajectory following this path would follow close to the line $y=1-x$ which corresponds to $\Lambda = 0$, this supports our choice of model as this trajectory relfects the path our universe is taking, in particular we observe a small and positive value for $\Lambda$.  

## Coordinate Transformations

We will now demonstrate a regular transformation applied to our dynamical system.  The aim here is to provide a concrete example of what a well-defined transformation should look like.

Our transformation will consist of two parts, first we will reflect along the line $y=x$ using the transformation $$V = \begin{pmatrix}
0& 1 \\
1 & 0
\end{pmatrix} $$and then we will apply a scaling in the $x$ and $y$ directions $$U = \begin{pmatrix}
x+2& 0 \\
0 & y+1
\end{pmatrix} $$We combine these two transformations via matrix multiplication and use the resulting matrix to define our transformed variables $$UV = \begin{pmatrix}
x+2& 0 \\
0 & y+1
\end{pmatrix}\begin{pmatrix}
0& 1 \\
1 & 0
\end{pmatrix} = \begin{pmatrix}
0& x+2 \\
y+1 & 0
\end{pmatrix} = T$$
$$ T\begin{pmatrix}
x \\
y
\end{pmatrix} = \begin{pmatrix}
m \\
r
\end{pmatrix}$$
Our critical points $R$, $O$, and $M$ get mapped to $(2,0)$, $(0,0)$, and $(0,1)$ respectively.  Our goal now is to apply this transformation to our dynamical system, we start with our expressions for $m$ and $r$ $$\begin{align*}
    m &= y(x+2)\\
    r &= x(y+1)
\end{align*}$$
Next we rearrange for $x$ and $y$ $$\begin{align*}
    x &= \frac{1}{2}(-2-m+r+\sqrt{8r+(-2-m+r)^2})\\
    y &= \frac{1}{2}(m-r+\frac{1}{2}(-2-m+r+\sqrt{8r+(-2-m+r)^2}))
\end{align*}$$
Now we can construct our system by differentiating $m$ and $r$ with respect to $t$ and sub in our equations for $x$, $x'$, $y$, and $y'$.

We finally can simulate this dynamical system on MatLab to create our new phase portrait

![alt text](https://github.com/BillGriffin98/CosmologicalDynamicSystems/blob/main/transformedtrianglepicture.png)

The determinant of our transformation of the Jacobian $J$ of our transformation $T$ is Det($J$) = $-(x+2)(y+1)$.  We evaulate this at each of our critical points $$\text{Det }J_R = -4,\phantom{....}\text{Det }J_O = -2,\phantom{....}\text{Det }J_M= -4 $$ hence since this is not singular at any point we expect all our critical points to be well-mapped to new critical points in the new dynamical system $(m',r')$.  We can see that our stable node at the origin $O$ has been unmoved by the mapping, while our saddle point $R$ and unstable node $M$ have been swapped and scaled by a factor of 2 in the $x$ direction.  These features present on our phase portrait reflect what we expected from our analysis.

## Irregular transformation
Now we will consider an alternate example that will demonstrate what can happen when our transformation is singular at certain points.\\

We will consider the following transformation:
$$T = \begin{pmatrix}
\frac{1}{1-y+\epsilon}& 0 \\
0 & \frac{1}{1-y+\epsilon}
\end{pmatrix}$$ 

which will scale our triangle by a factor of $\frac{1}{1+\epsilon}$ in the $x$ direction and a factor of $\frac{1}{\epsilon}$ in the $y$ direction.  Our critical points $R$, $O$, and $M$ are mapped to $(0,\frac{1}{\epsilon})$, $(0,0)$, and $(\frac{1}{\epsilon+1},0)$ respectively.

To construct our transformed dynamical system we proceed in a similar manor to before: $$\begin{pmatrix}
m\\r\end{pmatrix} = T\begin{pmatrix}
x\\y
\end{pmatrix}$$

We rearrange to get expressions for $x$ and $y$ in terms of $m$ and $r$ 
$$\begin{align*}
    x &= \frac{m(1+\epsilon)}{1+r}\\
    y &= \frac{r(1+\epsilon)}{1+r}
\end{align*}$$
Next we differentiate $m$ and $r$ and sub in our expressions for $x,x',y,y'$.$$\begin{align*}
    m' &= m(-3+3m(1+\epsilon)+4\epsilon r)\\
    r' &= r(-4+3m(1+\epsilon)+4\epsilon r)
\end{align*}$$
Now we can plot our phase portrait for various values of $\epsilon$.
\begin{figure}
\centering
\begin{subfigure}{.5\textwidth}
  \centering
  \includegraphics[width=1\linewidth]{epsilon1.png}
  \caption{$\epsilon=1$}
  \label{fig:sub1}
\end{subfigure}%
\begin{subfigure}{.5\textwidth}
  \centering
  \includegraphics[width=1\linewidth]{epsilon2.png}
  \caption{$\epsilon=1/2$}
  \label{fig:sub2}
\end{subfigure}
\label{fig:test}

\centering
\begin{subfigure}{.5\textwidth}
  \centering
  \includegraphics[width=1\linewidth]{epsilon10.png}
  \caption{$\epsilon=1/10$}
  \label{fig:test1}
\end{subfigure}%
\begin{subfigure}{.5\textwidth}
  \centering
  \includegraphics[width=1\linewidth]{epsilon100.png}
  \caption{$\epsilon=1/100$}
  \label{fig:test2}
 
\end{subfigure}
\caption{Varying values of $\epsilon$}
\end{figure}
We are interested in what happens to this system in the limit $\epsilon \rightarrow 0$.  When we evaluate the determinant of the Jacobian at each of our critical points we get $$\text{Det }J_R = (\frac{1}{\epsilon})^2,\phantom{....}\text{Det }J_O = (\frac{1}{1+\epsilon})^2,\phantom{....}\text{Det }J_M= (\frac{1}{1+\epsilon})^2 $$In the limit $\epsilon \rightarrow 0$ we can see Det $J_O$ = Det $J_M$ = 1 and hence are well-mapped, while Det $J_R$ tends to infinity is therefore singular.  We can see the issues with this transformation in our phase portraits where our triangle appears to be tending towards a rectangle.  With our standard dynamical system techniques we are unequipped to deal with a point that is mapped to infinity and hence the critical point $R$ is lost by the transformation.  We may have some hope of recovering it if we utilised techniques from projective geometry.  Hence we have demonstrated that a non-regular transformation will lose critical points been transforming between dynamical systems.

![alt text](https://github.com/BillGriffin98/CosmologicalDynamicSystems/blob/main/ChangesInEpsilon.png)

## $f(R)$ Gravity

The observation of a Type 1a supernova in 1998 provided sufficient evidence to deduce that the expansion of the universe is accelerating, this is one of many puzzling observations made in the last 100 years that have called for an extension to Einstein's theory of General Relativity (GR).  One such extension comes from generalising the Einstein-Hilbert action found at the foundation of GR $$\begin{align}
    S &= \frac{1}{2\kappa ^2}\int (f(R) + \mathcal{L}_m + \mathcal{L}_r)\sqrt{-\text{det } g}d^4x \\
\end{align}$$ 
where $\mathcal{L}_m, \mathcal{L}_r$ are the Lagrangian matter and radiation densities respectively.  For the remainder of this paper we will only be interested in vacuum solutions $(\mathcal{L}_m =  \mathcal{L}_r= 0 $) and hence our action reduces to $$\begin{align}
    S &= \frac{1}{2\kappa ^2}\int f(R)\sqrt{-\text{det } g}d^4x \\
\end{align}$$
As we can see the generalisation comes from substituting $R$ in the action for a more general $f(R)$, as a result this class of theories are called $f(R)$ theories.

We will be discussing two separate choices of variables and their resulting dynamical systems.  Both choices of variables are vacuum solutions of the FLRW universe and are derived from the same fundamental equations 
$$\begin{align}
H^2 &= \frac{R}{6} - \frac{f}{6F} - \frac{H\dot{F}}{F}\\
-2F\dot{H} &= \ddot{F} - \dot{F}H\\
R &= 6(2H^2 + \dot{H})
\end{align}$$

where $$F := \frac{\text{d}f}{\text{d}R}$$

## Amendola's method

Amendola begins by taking our field equation and dividing by $H^2$ $$1 =  \frac{R}{6H^2} - \frac{f}{6FH^2} - \frac{\dot{F}}{FH}.$$We can now introduce our dimensionless variables for our $f(R)$ model
$$\begin{align*} 
       x_{1} &= -\frac{\dot{F}}{FH} \\
       x_{2} &= -\frac{f}{6FH^2} \\
       x_{3} &= \frac{R}{6H^2} \\
\end{align*}$$
We can now rewrite our field equation in terms of these variables $$\begin{equation}
    1 = x_{1} + x_{2} + x_{3}\end{equation}$$
(this equation is from https://arxiv.org/pdf/0805.1726.pdf)
If we differentiate by $t$ we can get expressions for $\dot{x_1},...,\dot{x_4}$, we then convert into conformal time by applying the following transformation $$x' = \frac{1}{H}\dot{x}. $$We derive the following system: 
$$\begin{align}
    x_1' &= -1-x_3-3x_2+x_1^2-x_1 x_3 \\
    x_2' &=\frac{x_1 x_3}{m}-x_2(2x_3 - 4 - x_1) \\
    x_3' &= \frac{-x_1 x_3}{m}-2x_3(x_3-2)\\
    \end{align}$$
Where $m(r) = -r-1$ with $r = \frac{x_3}{x_2}$.  In our case of $f(R) = R + \alpha R^{-n}$ we can set $m = -n-1 = 1$.  We can apply to simplify our system down to 2 dimensions.
$$\begin{align}
    x_1' &= -4+2x_3+3x_1+x_1^2-x_1x_3 \\
    x_3' &= -x_1x_3^2 - 2x_3(x_3 - 2)
\end{align}$$
We can set these equations to 0 to obtain the following 3 critical points 

 $$\begin{equation} P_1: (0,2) \phantom{.....} P_2: (1,0),\phantom{.....} P_3: (-4,0)\end{equation}$$
 
 Now we can apply linear stability analysis to each critical point and assess the behaviour of the system near each point.  The Jacobian of the system is $$\begin{equation}
 \begin{pmatrix}
3+2x_1-x_3 & 2-x_1 \\
-x_3^2 & -2x_1x_3-4x_3+4
\end{pmatrix} 
 \end{equation}$$
 
From here we can find the find the Jacobian for each point
$$\begin{equation}
J_{P_1} = \begin{pmatrix}
1& 2 \\
-4 & -4
\end{pmatrix} \phantom{.....}J_{P_2} = \begin{pmatrix}
5& 1 \\
0 & 4
\end{pmatrix}\phantom{.....}J_{P_3} = \begin{pmatrix}
-5& 6 \\
0 &4
\end{pmatrix}\end{equation}$$

with trace and determinant $\{-3,4\}$, $\{-9,20\}$, and $\{-1,-20\}$ respectively.  Hence $P_1$ is a stable focus, $P_2$ is an unstable node, and $P_3$ is a saddle point.  We can plot the phase portrait of this system and verify this analytical analysis 

![alt text](https://github.com/BillGriffin98/CosmologicalDynamicSystems/blob/main/amendola.png)

## Ahlo's Method

Ahlo begins by attempting to simplify our field equation as must as possible.  First he substitutes in $f = R+\alpha R^2$ and $F = 1 + 2\alpha R$ $$\begin{equation}
    R^2 = 12H(\dot{R} + HR + \frac{H}{2\alpha})
\end{equation}$$
and then separates $$\begin{equation}
    R^2 = (\sqrt{\frac{12}{\alpha}}H)(\sqrt{12\alpha}(\dot{R} +HR + \frac{H}{2\alpha}))
\end{equation}$$
and then introduces the variables $$\begin{align}
    R^2 &= (T-X)(T+X)\\
    &= T^2 - X^2
\end{align}$$
which has greatly simplified our constraint.  Hence our choice of variables are
$$\begin{align}
                               H &= \sqrt{\frac{\alpha}{12}}(T-X) \\
    \dot{R}+HR+\frac{H}{2\alpha} & = \frac{1}{\sqrt{12\alpha}}(T+X) \\
\end{align}$$
with the constraint
$$\begin{equation} R^2 = T^2 - X^2
\end{equation}$$
                  
Note this constraint tells us that we are working in a 2D system.  From these equations we can solve for our dimensionless variables
$$\begin{align}
    T &= \sqrt{\frac{3}{\alpha}}(\alpha \dot{R} + \alpha HR + \frac{3}{2}H) \\
    X &= \sqrt{\frac{3}{\alpha}}(\alpha \dot{R} + \alpha HR - \frac{1}{2}H) \\
    R &= \sqrt{\frac{12}{\alpha}}(\alpha H\dot{R} + \alpha H^2R + \frac{1}{2}H^2)^{\frac{1}{2}} \\
\end{align}$$

Differentiating these variables gives us our dynamical system 
$$\begin{align*}
    \dot{T} &= \frac{1}{4\sqrt{3\alpha}}(R-2\alpha(T-X)^2)\\
    \dot{X} &= \frac{1}{4\sqrt{3\alpha}}(2\alpha(T-X)^2-3R)\\
    \dot{R} &= \frac{1}{4\sqrt{3\alpha}}(T+3X-2\alpha(T-X)R)
\end{align*}$$
Immediately we can apply (\ref{ahlo2d}) to obtain our 2D dynamical system, we take $\alpha = 1$ for simplicity.
$$\begin{align}
    \dot{T} &= \frac{1}{4\sqrt{3\alpha}}(\sqrt{T^2-X^2}-2(T-X)^2)\\
    \dot{X} &= \frac{1}{4\sqrt{3\alpha}}(2(T-X)^2-3\sqrt{T^2-X^2})
\end{align}$$
Next Ahlo introduces the following dimensionless variables $$\bar{X} = \frac{X}{T}, \phantom{.....}S= -\frac{R}{T}$$and hence our bound $R^2 = T^2 - X^2$ becomes $$\bar{X} + S = 1.$$ This constraint is globally satisfied via the following variables $$\bar{X} = \cos(\theta), \phantom{.....}S = \sin(\theta)$$ We also introduce the bounded variable
$$\begin{equation} \bar{T} = \frac{1}{1+2\alpha T}\end{equation}$$ 

giving us the bound $$0 < \bar{T} < 1.$$
Our equations now become the unconstrained system 
$$\begin{align*}
    \frac{\text{d}\bar{T}}{\text{d}\bar{t}} &= \bar{T}(1-\bar{T})[\bar{T}\sin(\theta) + (1-\bar{T})(1-\cos(\theta))^2] \\
    \frac{\text{d}\theta}{\text{d}\bar{t}} &= -\bar{T}(3+\cos(\theta))-(1-\bar{T})(1-\cos(\theta))\sin(\theta)
\end{align*}$$
where $\bar{t} = 2\sqrt{12\alpha}\bar{T}$.
This final system gives us two critical points of the form $(\bar{T},\theta)$ $$P_A: (0,2n\pi ), \phantom{.....} P_B: (0,\pi + 2n\pi)$$ 
and we can produce the Jacobian for this system 
$$\begin{equation*}
\begin{pmatrix}
(-1-\bar{T})(-1+3\bar{T})(-1+\cos(\theta))^2+\bar{T}\sin(\theta)(2-3\bar{T}) & \bar{T}(1-\bar{T})(\bar{T}\cos(\theta)+2\sin(\theta)(-1+\bar{T})(-1+\cos(\theta))) \\
-3+\sin(\theta)-\cos(\theta)(1+\sin(\theta)) & (-1+\bar{T})(\cos(\theta)-\cos(2\theta))+\bar{T}\sin(\theta)
\end{pmatrix} 
\end{equation*}$$
and hence we can find evaluate the Jacobian at each point
$$\begin{equation}
J_{P_A} = \begin{pmatrix}
0& 0 \\
-4 & 0
\end{pmatrix} \phantom{.....}J_{P_B} = \begin{pmatrix}
4& 0 \\
-2 & 2
\end{pmatrix}\end{equation}$$
Hence each point has the trace and determinant $\{0,4\}$ and $\{6,8\}$ respectively, corresponding to a saddle and an unstable node.


We can produce the phase portrait for this dynamical system 
![alt text](https://github.com/BillGriffin98/CosmologicalDynamicSystems/blob/main/Ahlo.png)

At this point, we can notice a discrepency in the number of critical points between Ahlo's and Amendola's dynamical systems. This is troubling as both systems were derived from the same fundemental equations, and hence should have an equal amount of critical points. 







