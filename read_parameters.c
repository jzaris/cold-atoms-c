#include "read_parameters.h"
#include "array2d.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

void initialize(const char *filename){
	FILE *f;
	char *line = NULL;
	size_t len = 0;

	f = fopen(filename, "r");  /*Input file with ensemble parameters*/
	if (f == NULL)
		exit(EXIT_FAILURE);

	int line_num = 0;
	int num_particles = 1;

	getline(&line, &len, f);
	int j = strlen(line);
	for (int i = 0; i<j; i++){
		if (line[i] == ' '){
			num_particles++;  /*Each space signifies another particle.  num_particles = 1 + num_spaces*/ 
		}
	}

	printf("Number of particles: %d \n", num_particles);
	double x[num_particles]; /*initial positions of each particle*/
	double y[num_particles];
	double z[num_particles];
	double vx[num_particles]; /*initial velocities*/
	double vy[num_particles];
	double vz[num_particles];
	double q[num_particles]; /*charge list*/
	double m[num_particles]; /*mass list*/
	double dt_arr[1]; /*timestep*/
	double duration_arr[1]; /*simulation duration*/
	double steps_arr[1]; /*number of timesteps*/
	double Bz_arr[1]; /*axial magnetic field*/
	double kz_arr[1]; /*stiffness parameter*/
	double delta_arr[1]; /*relative strength between rot. wall and trap potentials*/
        double omegaR_arr[1]; /*rotating wall frequency*/
        double phi0_arr[1]; /*initial rot. wall angle*/
        double dt;
	double duration;	
	double steps;
	double Bz;
	double kz;
	double delta;
	double omegaR;
	double phi0;
	double omegaB[num_particles];
	/*int steps;*/


	while (line_num<16){
		if (line_num==0){
			printf("%s", line);
			str_to_double(x, line, num_particles);  /*x-coordinates*/
			for (int i = 0; i < num_particles; i++)
				printf("%f ", x[i]);
			printf("\n");
		}

		if (line_num==1){
			getline(&line, &len, f); 
			printf("%s", line);
			str_to_double(y, line, num_particles);  /*y-coordinates*/
			for (int i = 0; i< num_particles; i++)
				printf("%f ", y[i]);
			printf("\n");
		} 

	        if (line_num==2){
                        getline(&line, &len, f);
                        printf("%s", line);
                        str_to_double(z, line, num_particles);  /*z-coordinates*/
                        for (int i = 0; i< num_particles; i++)
                                printf("%f ", z[i]);
                        printf("\n");
                }

		if (line_num==3){
			getline(&line, &len, f); 
			printf("%s", line);
			str_to_double(vx, line, num_particles);  /*x-velocities*/
			for(int i = 0; i<num_particles; i++)
				printf("%f ", vx[i]);
			printf("\n");
		}   
                if (line_num==4){
			getline(&line, &len, f);
			printf("%s", line);
			str_to_double(vy, line, num_particles);  /*y-velocities*/
			for(int i = 0; i<num_particles; i++)
			       printf("%f ", vy[i]);
			printf("\n");	
		}

		if (line_num==5){
                        getline(&line, &len, f);
                        printf("%s", line);
                        str_to_double(vz, line, num_particles);  /*z-velocities*/
                        for (int i = 0; i< num_particles; i++)
                                printf("%f ", vz[i]);
                        printf("\n");
                }
		if (line_num==6){
                        getline(&line, &len, f);
                        printf("%s", line);
                        str_to_double(q, line, num_particles);  /*charges*/
                        for(int i = 0; i<num_particles; i++)
                               printf("%f ", q[i]);
                        printf("\n");
                }

		if (line_num==7){
                        getline(&line, &len, f);
                        printf("%s", line);
                        str_to_double(m, line, num_particles);  /*masses*/
                        for(int i = 0; i<num_particles; i++)
                               printf("%f ", m[i]);
                        printf("\n");
                }


                if (line_num==8){
			 getline(&line, &len, f);  /*dt*/
			 printf("%s", line);
			 str_to_double(dt_arr, line, 1);
		         dt = dt_arr[0];
			 printf("%f ", dt);
			 printf("\n");
		}

		if (line_num==9){
                         getline(&line, &len, f);  /*duration*/
                         printf("%s", line);
                         str_to_double(duration_arr, line, 1);
			 duration = duration_arr[0];
                         printf("%f ", duration);
			 printf("\n");
                }

		if (line_num==10){
			 getline(&line, &len, f);  /*number of steps*/
                         printf("%s", line);
			 str_to_double(steps_arr, line, 1);
			 steps = steps_arr[0];
			 printf("%f ", steps);
                         printf("\n");	 
		}

		if (line_num==11){
                         getline(&line, &len, f);  /*axial magnetic field*/
                         printf("%s", line);
                         str_to_double(Bz_arr, line, 1);
                         Bz = Bz_arr[0];
                         printf("%f ", Bz);
                         printf("\n");
                }

                if (line_num==12){
                         getline(&line, &len, f);  /*kz*/
                         printf("%s", line);
                         str_to_double(kz_arr, line, 1);
                         kz = kz_arr[0];
                         printf("%f ", kz);
                         printf("\n");
                }

		if (line_num==13){
                         getline(&line, &len, f);  /*delta*/
                         printf("%s", line);
                         str_to_double(delta_arr, line, 1);
                         delta = delta_arr[0];
                         printf("%f ", delta);
                         printf("\n");
                }

		if (line_num==14){
                         getline(&line, &len, f);  /*omegaR*/
                         printf("%s", line);
                         str_to_double(omegaR_arr, line, 1);
                         omegaR = omegaR_arr[0];
                         printf("%f ", omegaR);
                         printf("\n");
                }

		if (line_num==15){
                         getline(&line, &len, f);  /*phi0*/
                         printf("%s", line);
                         str_to_double(phi0_arr, line, 1);
                         phi0 = phi0_arr[0];
                         printf("%f ", phi0);
                         printf("\n");
                }

		line_num = line_num+1;
	}
	

	arr.Initialize(9, num_particles); /*rows correspond to pos (x,y,z), vel (x,y,z), q, m, omegaB*/
	for(int i = 0; i<num_particles; i++){
                arr(0,i) = x[i];
        }
        for(int i = 0; i<num_particles; i++){
                arr(1,i) = y[i];
        }
	for(int i = 0; i<num_particles; i++){
                arr(2,i) = z[i];
        }
        for(int i = 0; i<num_particles; i++){
                arr(3,i) = vx[i];
        }
        for(int i = 0; i<num_particles; i++){
                arr(4,i) = vy[i];
        }
	for(int i = 0; i<num_particles; i++){
                arr(5,i) = vz[i];
        }
        for(int i = 0; i<num_particles; i++){
                arr(6,i) = q[i];
        }
        for(int i = 0; i<num_particles; i++){
                arr(7,i) = m[i];
        }

	/*omega_B for each particle*/
	for(int i = 0; i<num_particles; i++){
		omegaB[i] = Bz *q[i] /m[i];
                arr(8,i) = omegaB[i];
        }

	for (int i = 0; i < num_particles; i++)
        	printf("%f ", omegaB[i]);
        printf("\n");

	/*assign values to global variables*/
	dt_global = dt;
        duration_global = duration;
        steps_global = steps;
        particles_global = num_particles;
        Bz_global = Bz;
        kz_global = kz;
        delta_global = delta;
        omegaR_global = omegaR;
        phi0_global = phi0;
        phi_global = phi0;


	//initialize forces array with all zeros
	forces.Initialize(3, num_particles);
	for(int i = 0; i<3; i++){
        	for(int j = 0; j<num_particles; j++){
                	forces(i,j) = 0.0;
        	}
	}
}



void str_to_double(double *arr, const char *str, int num_particles){   /*takes a string (str) of (num_particles) doubles and fills an array (arr) with the doubles*/
	int math_counter = 0;  /*indexes the digits in the mantissa*/
	int exp_counter =0;    /*indexes the digits in the exponent*/
	int array_counter =0;  /*indexes the number of particles*/
	double this_val = 0;   /*holds value of mantissa*/
	int this_exp = 0;      /*holds value of exponent*/
	bool pre_exp = true;   /*true if the current digit is part of the mantissa, false if it is part of the exponent*/
	bool exp_neg = false;  /*true if the exponent is negative*/
	bool val_neg = false;  /*true if the mantissa is negative*/
	char this_char;
	for(int i = 0; i < strlen(str); i++){   /*indexes each character*/
		/*printf("%c ", str[i]);*/
		/*printf("%f ", this_val);*/
		if(str[i] == ' '){             /*space delimits the end of a value*/
			/*printf("%f ", this_val);
			printf("%i ", this_exp);*/
			if(val_neg){          
				if(exp_neg)
					arr[array_counter] = -1*this_val*pow(10,-1*this_exp);  /*add value, with mantissa 'this_val' and exponent 'this_exp' to arr*/
				else
					arr[array_counter] = -1*this_val*pow(10, this_exp);
			}
			else{
				if(exp_neg)
					arr[array_counter] = this_val*pow(10, -1*this_exp);
				else
					arr[array_counter] = this_val*pow(10, this_exp);
			}

			array_counter++;   /*reset to default values before reading next value*/
			math_counter = 0;
			exp_counter =0;
			this_val = 0;
			this_exp = 0;
			exp_neg = false;
			val_neg = false;
			pre_exp = true;
		}

		else{
			if(pre_exp)
				pre_exp = !(str[i] == 'e');    /*check to see if current character is e*/
			
			if (pre_exp){
				if(math_counter == 0){     
					if(str[i] == '-'){
						val_neg = true;   /*if the first character of the value is the minus sign, then the value is negative*/
					}
					else{
						this_char = str[i];
						this_val += atoi(&this_char);  /*the first numeric character occupies the ones place of the mantissa*/
						math_counter++;
						/*printf("%c ", str[i]);
						printf("%i ",atoi(&str[i]));
						printf("%f ", this_val);*/
					}	
					
				}
				else{
					if(math_counter ==1)        /*this character is the decimal point*/
						math_counter++;
		
					else{
						this_char = str[i];
						this_val += atoi(&this_char)*pow(10, -(math_counter-1));   /*each digit of the mantissa after the decimal point*/
						math_counter++;
					}
				}
			}
			else{  /*means current index corresponds to the exponent*/

				if(str[i] == '-' || str[i] == '+'){
					if (str[i] == '-')
						exp_neg = true;     /*if the first character of the exponent is a minus sign, then the exponent is negative*/
				}

				else{
					if(str[i] != 'e'){  /*ignore e, it doesn't provide any info*/
						if(exp_counter == 0){
							this_char = str[i];
							this_exp += 10*atoi(&this_char); /*the exponent has two digits, the first of which is the tens digit*/
							exp_counter++;		
						}
						else{
							this_char = str[i];
							this_exp += atoi(&this_char);  /*ones digit*/
						}
						/*printf("This exp %i ", this_exp);*/
					}

				}
					
			}

		}
		
		if(i == strlen(str) - 1){    /*denotes the end of the final value*/
			/*printf("%f ", this_val);
			printf("%i ", this_exp);*/
			if(val_neg){
				if(exp_neg)
					arr[array_counter] = -1*this_val*pow(10,-1*this_exp);   /*add final value to arr*/
				else
					arr[array_counter] = -1*this_val*pow(10, this_exp); 
			}	
			else{
				if(exp_neg)
					arr[array_counter] = this_val*pow(10, -1*this_exp);
				else
					arr[array_counter] = this_val*pow(10, this_exp); 
			}
		}
	}
	

}
