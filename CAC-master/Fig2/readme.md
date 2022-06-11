Please run main.m to get the CAC landscape in Figure 2.
Force5i.m records the CAC model (ODEs), 
multivariate_normal_distribution is the density function of the Gaussian distribution which is used to calculate the density function of CAC model;
Solver.m is used to solve the ODEs, calculate the mean value and the covariance of each stable state and calculate the transition paths and actions between the stable states.
jacobian.m is the jacobian matrix of the CAC model, 
plot_Landscape.m is the main function which is uesd to calculate the density function of expression level of the system and use the DRL method [1] to plot the landscape and transition paths of the system.
Please run plot_Landscape.m and the program runs about 50 minutes to get landscape, tranistion path and stationary distribution of stable states.