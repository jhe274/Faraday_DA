import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from constants import Constants as Consts
from theory import Theory

class Manipulate:

    def __init__(self):
        self.consts = Consts()
        self.theory = Theory()
    
    def manipulate_plot(self, Kn, T, B, P, gamma_D1, gamma_D2, const):
        nu = np.linspace(self.consts.Nu39_D2 - 10 * 1e9, self.consts.Nu39_D2 + 10 * 1e9, 8000)
        y_min0, y_max0 = -200, 200

        # Create the figure and the line that we will manipulate
        fig, ax = plt.subplots()
        plt.subplots_adjust(left=.1, bottom=.4)
        y = self.theory.resonant_FR(nu, Kn, T, B, P, gamma_D2, gamma_D2, const) * 1e6                                       # [μrad]
        l, = plt.plot(nu * 1e-9 - self.consts.Nu39_D2 * 1e-9, y, lw=2)
        ax.set_xlabel('Frequency (GHz)')
        ax.set_ylabel('Faraday Rotation (μrad.)')
        ax.set_ylim(y_min0, y_max0)
        
        # Define sliders for each parameter
        ax_a = plt.axes([.1, .3, .35, .03])
        ax_b = plt.axes([.1, .25, .35, .03])
        ax_c = plt.axes([.1, .2, .35, .03])
        ax_d = plt.axes([.1, .15, .35, .03])
        ax_e = plt.axes([.1, .1, .35, .03])
        ax_f = plt.axes([.1, .05, .35, .03])
        ax_g = plt.axes([.55, .3, .35, .03])
        ax_ymin = plt.axes([.55, .25, .35, .03])
        ax_ymax = plt.axes([.55, .2, .35, .03])

        slider_a = Slider(ax_a, 'Kn', .1, 10, valinit=Kn)
        slider_b = Slider(ax_b, 'T', 10, 200, valinit=T)
        slider_c = Slider(ax_c, 'B', -40, 40, valinit=B)
        slider_d = Slider(ax_d, 'P', -1, 1, valinit=P)
        slider_e = Slider(ax_e, '$\gamma_{D1}$', 0, 2000, valinit=gamma_D1)
        slider_f = Slider(ax_f, '$\gamma_{D2}$', 0, 2000, valinit=gamma_D2)
        slider_g = Slider(ax_g, 'Offset', -500, 500, valinit=const)
        slider_ymin = Slider(ax_ymin, 'y_min', -5000, 0, valinit=y_min0)
        slider_ymax = Slider(ax_ymax, 'y_max', 0, 5000, valinit=y_max0)

        # Update function to modify the plot
        def update(val):
            a = slider_a.val
            b = slider_b.val
            c = slider_c.val
            d = slider_d.val
            e = slider_e.val
            f = slider_f.val
            g = slider_g.val
            y_min = slider_ymin.val
            y_max = slider_ymax.val

            y = self.FR(nu, a, b, c, d, e, f, g)
            l.set_ydata(y)
            ax.set_ylim(y_min, y_max)
            fig.canvas.draw_idle()

        # Connect sliders to update function
        slider_a.on_changed(update)
        slider_b.on_changed(update)
        slider_c.on_changed(update)
        slider_d.on_changed(update)
        slider_e.on_changed(update)
        slider_f.on_changed(update)
        slider_g.on_changed(update)
        plt.show()

if __name__ == "__main__":
     plotter = Manipulate()

     plotter.manipulate_plot(1.5, 21.7, -5, 0, 6, 6, 50)