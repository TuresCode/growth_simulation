import streamlit as st
import numpy as np
import matplotlib.pyplot as plt

# Define the default parameter values

params = {
    'k1': (0.7, 'Maximum growth rate on glucose (h^-1)'),
    'k2': (0.3, 'Maximum growth rate on lactose (h^-1)'),
    'Ks1': (0.1, 'Glucose concentration at which the growth rate is half its maximum (g/L)'),
    'Ks2': (0.2, 'EThanol concentration at which the growth rate is half its maximum (g/L)'),
    'Y1': (0.5, 'Yield coefficient for glucose (g cells/g glucose)'),
    'Y2': (0.3, 'Yield coefficient for ethanol (g cells/g lactose)'),
    'S0': (12, 'Initial substrate concentration (g/L)'),
    'V0': (0.001, 'Initial cell concentration (g/L)'),
    'glucose_to_ethanol': (1, 'Conversion factor for glucose to ethanol'),
    'lag_time': (50, 'Lag phase from glucose to ethanol')
}



# Create a Streamlit sidebar to allow the user to adjust the parameter values
st.sidebar.write("# Parameters")
params['k1'] = st.sidebar.slider('k1', min_value=0.0, max_value=1.0, value=params['k1'][0], step=0.05)
params['k2'] = st.sidebar.slider('k2', min_value=0.0, max_value=1.0, value=params['k2'][0], step=0.05)
params['Ks1'] = st.sidebar.slider('Ks1', min_value=0.0, max_value=1.0, value=params['Ks1'][0], step=0.05)
params['Ks2'] = st.sidebar.slider('Ks2', min_value=0.0, max_value=1.0, value=params['Ks2'][0], step=0.05)
params['Y1'] = st.sidebar.slider('Y1', min_value=0.0, max_value=1.0, value=params['Y1'][0], step=0.05)
params['Y2'] = st.sidebar.slider('Y2', min_value=0.0, max_value=1.0, value=params['Y2'][0], step=0.05)
params['S0'] = st.sidebar.slider('S0', min_value=0, max_value=50, value=params['S0'][0], step=1)
params['V0'] = st.sidebar.slider('V0', min_value=0.0, max_value=1.0, value=params['V0'][0], step=0.005)
params['glucose_to_ethanol'] = st.sidebar.slider('glucose_to_ethanol', min_value=0, max_value=10, value=params['glucose_to_ethanol'][0], step=1)
params['lag_time'] = st.sidebar.slider('lag_time', min_value=0, max_value=50, value=params['lag_time'][0], step=1)

st.write("# Diauxic Growth Simulation of Yeast on Glucose")
st.markdown('*S. cerevisiae*, as a Crabtree-positive yeast, predominantly ferments pyruvate to ethanol in high glucose conditions. When glucose or other preferred carbon sources are depleted, *S. cerevisiae* switches to aerobic respiration and utilizes ethanol as carbon source instead, a phenomenon known as diauxic shift.')



# Run the simulation
# Define the simulation function
@st.cache_resource
def growth_simulation(params):
    # Unpack the parameters
    k1 = params['k1']
    k2 = params['k2']
    Ks1 = params['Ks1']
    Ks2 = params['Ks2']
    Y1 = params['Y1']
    Y2 = params['Y2']
    S0 = params['S0']
    V0 = params['V0']
    glucose_to_ethanol = params['glucose_to_ethanol']
    lag_time = params['lag_time']
    end_time = 'Not calculated'

    # Time vector
    t = np.linspace(0, 20, 1000) # time in min
    # Glucose concentration
    S1 = np.zeros_like(t)
    S1[0] = S0
    # Lactose concentration
    S2 = np.zeros_like(t)

    # Simulation
    V = np.zeros_like(t)
    V[0] = V0
    substrate = "glucose"
    curr_time = 0
    check_end = True
    for i in range(1, len(t)):
        if S1[i-1] > 0.01:
            mu = k1*S1[i-1]/(Ks1 + S1[i-1])
            dVdt = mu*V[i-1] #* (mu2*V[i-1])  #growth rate * current cell concentration
            dS1dt = -mu/Y1*V[i-1] #yiel coefficient (how many new cells in gram per gram glucose)
            #dS2dt = -mu/Y2*V[i-1]
            V[i] = V[i-1] + dVdt*(t[i]-t[i-1]) #old cell concentration + growth rate * time interval
            S1[i] = S1[i-1] + dS1dt*(t[i]-t[i-1]) #how much substrate depleted
            S2[i] = (S0 - S1[i]) * glucose_to_ethanol
        
        elif curr_time <= lag_time:
            curr_time += 1
            V[i] = V[i - 1]
            S1[i] = S1[i - 1]
            S2[i] = S2[i - 1]
        else:
            mu = k2*S2[i-1]/(Ks2 + S2[i-1])
            dVdt = mu*V[i-1] #* (mu2*V[i-1])  #growth rate * current cell concentration
            dS2dt = -mu/Y2*V[i-1]
            V[i] = V[i-1] + dVdt*(t[i]-t[i-1]) #old cell concentration + growth rate * time interval
            S1[i] = S1[i-1] #S1[i-1] + dS1dt*(t[i]-t[i-1]) #how much substrate depleted
            S2[i] = S2[i-1] + dS2dt*(t[i]-t[i-1])
            
            if S2[i] < 0.01:
                V[i] = V[i-1]
                S1[i] = S1[i - 1]
                S2[i] = S2[i - 1]
                if check_end:
                    end_time = str(round(i/60,2))
                    check_end = False

      #print(V.argmax())
 

    return V, S1, S2, end_time

#if st.button('Run Simulation:'):
#    with st.spinner('Running my simulation...'):

        # Run the simulation
V, S1, S2, end_time = growth_simulation(params)
# plot the results
fig = plt.figure(figsize=(7, 7))
plt.plot(V, label="cells")
plt.plot(S1, label="glucose")
plt.plot(S2, label="ethanol")
plt.legend()
plt.xlabel("Time (min)")
#st.success('Done!')
st.pyplot(fig)

# Display the parameter values
st.write("## Current Parameter Values:")
st.write(params)
st.write('Time until final OD600 in hours: '+end_time)
st.write('')
st.write('Provided by Fabian & David 😀')
