import pandas as pd
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
import tensorflow as tf

# Load data into a Pandas DataFrame
df = pd.read_csv('demo/covariables_ml/data.csv')
print(df)
# csv file:
# paciente,V_vc,Cl_vc,Cph_vc,V_nac,Cl_nac,V_dfo,Cl_dfo,tac,sexo,edad,masa_corporal,talla,imc,tbq,oh,unidades_oh_semana
# p01,12.02860795,10.57049508,0.095697327,11.83961152,14.91912234,51.30713959,176.9365092,tac1,f,28,56,165,20.57,0,1,2
# p03,27.22422204,24.85702137,0.002027429,9.481120625,13.31978504,79.90111182,231.1407025,tac1,f,29,63,170,21.8,0,1,1
# p05,8.215491511,9.783662146,0.097026283,11.37407603,11.89733771,72.14017292,224.5281729,tac1,m,25,82.3,187,23.54,1,1,6
# p06,7.309362509,10.98943634,0.119753068,12.91065353,25.01807896,130.6309568,452.8731886,tac1,m,26,65,172,21.97,0,1,2
# p08,9.138905604,6.575794592,0.03866955,10.45385004,7.975516473,112.3993323,247.4420741,tac1,f,21,63,157,25.56,0,1,2
# p09,13.08837342,14.2021287,0.051329082,19.36253977,21.7196139,126.793044,211.4559219,tac1,m,31,96,176,30.99,1,1,5
# p10,7.768215334,7.35114482,0.133260311,11.47316366,22.47986455,151.8998961,461.3372656,tac2,f,27,61,162,23.24,0,1,1
# p11,14.74593003,7.854471222,0.017848161,32.0487653,16.99118496,67.82522704,223.8048188,tac2,m,23,62,165,22.77,0,0,0
# p14,13.56561988,11.59093048,0.078726288,5.445222391,28.7222153,126.8889788,544.7653052,tac2,m,23,87,191,23.85,0,0,0
# p15,8.759511897,7.560556239,0.014562463,19.54088788,22.20254161,398.535391,574.192524,tac2,m,21,61,164,22.68,0,0,0
# p17,13.68425645,15.55176645,0.082575209,7.547631289,73.52917886,371.5777454,856.3343094,tac2,m,27,98,183,29.26,0,1,0.25
# p18,12.87326526,10.89319418,0.049893316,10.43573105,11.66648512,247.6927164,734.7437109,tac2,f,23,55,171,18.81,0,1,3


# analyze how column V_vc from df is affected by the rest of the column from df with multiple linear regression
# define the independent variables
X = df[['Cl_vc', 'Cph_vc', 'V_nac', 'Cl_nac', 'V_dfo', 'Cl_dfo', 'edad', 'masa_corporal', 'talla', 'imc', 'tbq', 'oh',
        'unidades_oh_semana']]
# define the dependent variable
y = df['V_vc']
# create a linear regression model
model = LinearRegression()
# fit the model to the data
model.fit(X, y)
# get the coefficient of determination (R^2) of the prediction
r_sq = model.score(X, y)
print('coefficient of determination:', r_sq)
# get the intercept of the prediction
print('intercept:', model.intercept_)
# get the slope of the prediction
print('slope:', model.coef_)
# plot the data and the prediction
plt.scatter(X['Cl_vc'], y, color='black')
plt.scatter(X['Cph_vc'], y, color='red')
plt.scatter(X['V_nac'], y, color='green')
plt.scatter(X['Cl_nac'], y, color='blue')
plt.scatter(X['V_dfo'], y, color='yellow')
plt.scatter(X['Cl_dfo'], y, color='orange')
plt.scatter(X['edad'], y, color='purple')
plt.scatter(X['masa_corporal'], y, color='brown')
plt.scatter(X['talla'], y, color='pink')
plt.scatter(X['imc'], y, color='gray')
plt.scatter(X['tbq'], y, color='cyan')
plt.scatter(X['oh'], y, color='magenta')
plt.scatter(X['unidades_oh_semana'], y, color='olive')
plt.plot(X, model.predict(X), color='blue', linewidth=3)
plt.show()

