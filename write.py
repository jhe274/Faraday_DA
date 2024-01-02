import os
import analyze

class Write:

    def Theta_and_lambda_data_to_txt(self, date, lambda_path, lock_in_path, i, n):
        '''
            Append measured polarization rotation and wavelength data to csv file
        '''
        file_name = "Polarization_rotations_vs_wavelength.csv"
        N, A_fit, Lambda_ave, Lambda_std, Lambda_ste, Theta_ave, Theta_std, Theta_ste = analyze.filtered_theta_and_lambda(lambda_path, lock_in_path, i, n)

        # Create a list of data points
        data = [date, N, A_fit, Lambda_ave, Lambda_std, Lambda_ste, Theta_ave, Theta_std, Theta_ste]
        
        self.save_data_to_file(file_name, data)

    def save_data_to_file(self, path, filename, data):
        try:
            file_path = os.path.join(path, filename)
            header = "Date (MM-DD-YYYY), Number of Data Points, Fitted_scan_amplitude (m), lambda_ave (m), lambda_std (m), lambda_ste (m), theta_ave (rad), theta_std (rad), theta_ste (rad)\n"
            if not os.path.isfile(file_path):
                with open(file_path, "w") as file:
                    file.write(header)

            # Check for duplicate data
            with open(file_path, "a+") as file:
                file.seek(0)
                existing_data = file.readlines()
                for line in existing_data:
                    existing_columns = line.strip().split(",")
                    new_columns = [str(item) for item in data]
                    if existing_columns == new_columns:
                        print("Duplicate data detected, abort writing.")
                        return
            
            # Append new data
            with open(file_path, "a") as file:
                file.write(",".join(map(str, data)) + "\n")

            print("Data appended to the file successfully.")
        except Exception as e:
            print(f"An error occurred while saving data to the file: {e}")