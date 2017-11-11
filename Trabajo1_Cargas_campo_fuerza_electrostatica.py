
# coding: utf-8

# # Cargas, campo y fuerza electrostática
# ## Primer ejercicio
# Un electrón con una capacidad de trabajo de $72090×10^{−19}$ J orbita de manera perpendicular a un campo magnético de 3250 G. ¿Cuál es el radio de la órbita? ¿Y su frecuencia y periodo angular? Resuelve el mismo ejercicio para un antineutrón.
# ![](Alphaparticlemagnetic.png)
# 
# 

# ### Importar frameworks

# In[94]:


from pint import UnitRegistry
pintunits = UnitRegistry()
import sympy as sp
from sympy import *
import sympy.physics.units as sp_u
import scipy.constants as sp_c
import math

from sympy.physics.units import *
from sympy import *
import sympy.physics.units as pu


# ### Convertir de Gauss a Tesla

# In[95]:


valor = (3250 * pintunits.gauss).to(pintunits.tesla).magnitude * sp_u.tesla
valor


# ### Ecuaciones
# #### Energía cinética
# $E_c = \frac{1}{2}mv^2$
# #### Lorentz
# $q_e v_\mathrm{e}B = \frac{m_ev_e^2}{r}$

# In[96]:


# Declaramos las variables de todas las ecuaciones
q_e, v, v_e, B, m, m_e, r, Ec = symbols('q_e v v_e B m m_e r Ec', positive = True, real = True)
eq_lorentz = Eq(q_e*v_e*B, (m_e*v_e**2)/r)
e_cinetica = Eq(Ec, (m*v**2)/2)

ecuacion = (eq_lorentz.subs(v_e, solve(e_cinetica, v)[0])).subs(m, m_e)
ecuacion


# #### Resolución de la ecuación
# $\frac{\sqrt{2E_c}Bq_e}{\sqrt{m_e}} = \frac{2E_c}{r}$

# In[97]:


expresion = solve(ecuacion, r)[0]
expresion


# $r = \frac{\sqrt{2E_cm}}{Bq}$

# In[98]:


valor = 0.325 * sp_u.tesla
energia_cinetica = 72090E-19 * sp_u.joule
masa_electron = sp_c.electron_mass * sp_u.kg
carga_electron = sp_c.electron_volt * sp_u.coulombs
radio = expresion.subs([(B, valor), (m_e, masa_electron), (q_e, carga_electron), (Ec, energia_cinetica)])


# ### Radio de la órbita:

# In[99]:


N(radio)


# ## Frecuencia y Periodo angular

# $W = \frac{2pi}{T}$

# In[100]:


velocidad = solve(e_cinetica, v)[0]
velocidad


# $W = \frac{2pi}{2pir}{v}$

# In[101]:


tiempo = 2*pi*radio/velocidad
N(tiempo)


# In[102]:


frecuencia = 2*pi/tiempo
N(frecuencia)


# In[103]:


frecuencia2 = velocidad/radio
N(frecuencia2)


# ###### Antineutrón: no orbita. carga=0

# ## Segundo ejercicio
# Calcula el módulo de la fuerza magnética que actúa sobre un electrón proveniente del Sol que penetra en la aurora boreal joviana. Haz cálculos aproximados basados en la búsqueda de información relativa a Júpiter, su campo magnético y el fundamento físico de una aurora boreal. Asume que la velocidad del electrón es prácticamente la de la luz.
# ![](http://en.es-static.us/upl/2011/08/jupiter-aurora.jpg)
# Fuente: http://en.es-static.us/upl/2011/08/jupiter-aurora.jpg

# #### Lorentz force
# 
# $$F = {qv}{B}$$

# q: carga $q = 1.6*10^{-19}$ C <br />
# v: velocidad de partícula (velocidad de la luz) $v = 3*10^{8}$ m/s <br />
# B: campo magnético $B = 778*10^{6}$ Km 

# In[104]:


F, q, v, B = symbols('F q v B', positive = True, real = True)
f_lorentz = Eq(F, (q*v*B))
ecu = solve(f_lorentz, F)[0]
ecu


# In[105]:


carga = 1.6E-19 * sp_u.coulombs
velocidad_particula = sp_c.speed_of_light * sp_u.m / sp_u.second
campo_magnetico = 778E6 * sp_u.km

fuerza_magnetica = ecu.subs([(q, carga), (v, velocidad_particula), (B, campo_magnetico)])
N(fuerza_magnetica)


# ## Tercer ejercicio
# 
# Se tiene una cantidad pequeña de material lubricante de masa 2.41×1010 u (unidades de masa atómica) y una carga de 4.8×10−19 C. La gota de aceite se encuentra flotando en equilibrio gracias a la harmonía de la fuerza gravitatoria más otra fuerza extra de naturaleza eléctrica. ¿Cuál es la dirección y magnitud del campo eléctrico originado por dicha fuerza? ¿A qué te recuerda la experiencia descrita? Justifica tu respuesta.
# ![](t2actividad3.png)

# $$E = \frac{mg}{q}$$

# $m = 2.41*10^{10}$ U <br />
# $g = 9.81$ N <br />
# $q = 4.8*10^{-19}$ C 

# In[106]:


m, g, q, E = symbols('m g q E', positive = True, real = True)

m = 2.41E10 * sp_u.amu.convert_to(pu.kilogram)
g = sp_u.acceleration_due_to_gravity
q = 4.8E-19 * sp_u.C

formula = Eq(E, ((m*g)/q))
resultado = convert_to((solve(formula, E)[0]), volt/meter)
N(resultado)

