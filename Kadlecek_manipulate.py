import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from constants import Constants as Consts
from theory import Theory

class Plot:

    def __init__(self):
        self.consts = Consts()
        self.theory = Theory()

    def FR(self, nu, Kn, T, B, P):
            l = (7.5-0.159 * 2) * 1e-2 
            delta_nu_D2 = nu - self.consts.Nu39_D2
            delta_nu_D1 = nu - self.consts.Nu39_D1
            delta_doppler_D2 = self.theory.doppler_broad(self.consts.Nu39_D2, T)
            delta_doppler_D1 = self.theory.doppler_broad(self.consts.Nu39_D1, T)

            dia_FR1 = ( 
                (7*(delta_nu_D2**2 - delta_doppler_D2**2/4) / (delta_nu_D2**2 + delta_doppler_D2**2/4)**2) + 
                (4*(delta_nu_D1**2 - delta_doppler_D1**2/4) / (delta_nu_D1**2 + delta_doppler_D1**2/4)**2) - 
                (2 / (delta_nu_D2 * delta_nu_D1))) / (3 * self.consts.h)
                # (2*(delta_nu_D1*delta_nu_D2) / ((delta_nu_D2-delta_doppler_D2) * (delta_nu_D1-delta_doppler_D1))**2)) / (3 * self.consts.h)
            
            # dia_FR2 = np.sign(B) * ((nu / (self.consts.Nu39_D1 * (nu - self.consts.Nu39_D1))) 
            #                            - (nu / (self.consts.Nu39_D2 * (nu - self.consts.Nu39_D2)))) / (self.consts.k_B * (273.15 + T))
            
            dia_FR = self.consts.mu_B * np.sign(B) * B * 1e-4 * (dia_FR1)

            para_FR = P * (
                    (delta_nu_D2 / ((delta_nu_D2 - delta_doppler_D2) ** 2 + self.theory.doppler_broad(self.consts.Nu39_D2, T)**2/4)) -
                    (delta_nu_D1 / ((delta_nu_D1 - delta_doppler_D1) ** 2 + self.theory.doppler_broad(self.consts.Nu39_D1, T)**2/4))
                )
            theta = self.consts.alpha * Kn * 1e14 * l * (dia_FR + para_FR) * 1e6                                                  # [microrad]
            return theta
    
    def manipulate_plot(self, Kn, T, B, P):
        nu = np.linspace(self.consts.Nu39_D2 - 10 * 1e9, self.consts.Nu39_D2 + 10 * 1e9, 8000)
        y_min0, y_max0 = -200, 50

        # Create the figure and the line that we will manipulate
        fig, ax = plt.subplots()
        plt.subplots_adjust(left=0.1, bottom=0.4)
        y = self.FR(nu, Kn, T, B, P)
        l, = plt.plot(nu * 1e-9 - self.consts.Nu39_D2 * 1e-9, y, lw=2)
        ax.set_xlabel('Frequency (GHz)')
        ax.set_ylabel('Faraday Rotation (microrad.)')
        ax.set_ylim(y_min0, y_max0)
        

        # Define sliders for each parameter
        ax_a = plt.axes([0.25, 0.3, 0.65, 0.03])
        ax_b = plt.axes([0.25, 0.25, 0.65, 0.03])
        ax_c = plt.axes([0.25, 0.2, 0.65, 0.03])
        ax_d = plt.axes([0.25, 0.15, 0.65, 0.03])
        ax_ymin = plt.axes([0.25, 0.1, 0.65, 0.03])
        ax_ymax = plt.axes([0.25, 0.05, 0.65, 0.03])

        slider_a = Slider(ax_a, 'Kn', .1, 10, valinit=Kn)
        slider_b = Slider(ax_b, 'T', 20, 23, valinit=T)
        slider_c = Slider(ax_c, 'B', -40, 40, valinit=B)
        slider_d = Slider(ax_d, 'P', -1, 1, valinit=P)
        slider_ymin = Slider(ax_ymin, 'y_min', -5000, 0, valinit=y_min0)
        slider_ymax = Slider(ax_ymax, 'y_max', 0, 5000, valinit=y_max0)

        # Update function to modify the plot
        def update(val):
            a = slider_a.val
            b = slider_b.val
            c = slider_c.val
            d = slider_d.val
            y_min = slider_ymin.val
            y_max = slider_ymax.val

            y = self.FR(nu, a, b, c, d)
            l.set_ydata(y)
            ax.set_ylim(y_min, y_max)
            fig.canvas.draw_idle()

        # Connect sliders to update function
        slider_a.on_changed(update)
        slider_b.on_changed(update)
        slider_c.on_changed(update)
        slider_d.on_changed(update)

        plt.show()

if __name__ == "__main__":
     plotter = Plot()

     plotter.manipulate_plot(1.5, 21.7, -5, 0)