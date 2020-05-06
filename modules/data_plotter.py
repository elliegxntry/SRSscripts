import numpy as np
import os
import data_wrangler as dw
from matplotlib.patches import Circle
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import AxesGrid
from matplotlib.ticker import LogFormatterSciNotation
from matplotlib.ticker import FixedLocator
import h5py

class data_plotter(dw.data_wrangler):

    def __init__(self, config, project, specs, to_load=False, **kwargs):
        dw.data_wrangler.__init__(self, config, project, specs, to_load, **kwargs)
        plt.rc('font', size=12)
        plt.rc('axes', labelsize=16)
        plt.rc('xtick', labelsize=12)
        plt.rc('ytick', labelsize=12)
        plt.rc('legend', fontsize=12)
        plt.rc('figure', titlesize=12)

    #------------------------------------------------------------------------
    def plot_vavg_hst(self, time, qs=None, **kwargs):
        """
        options:
        labels,
        legend labels (leg_lab)
        ylim
        xlim
        title
        linestyles
        colors
        """
        ls = kwargs.get("ls", [":", "--", "-", "-."])
        xlim = kwargs.get("xlim", [time[0], time[-1]])
        ylim = kwargs.get("ylim", [])
        xlabel = kwargs.get("xlabel", "Time (orbits)")
        ylabel = kwargs.get("ylabel", None)
        colors = kwargs.get("colors")
        leg_lab = kwargs.get("leg_lab", list(qs.keys()))
        norm = kwargs.get("norm", True)
        ylog = kwargs.get("ylog", False)
        normval = kwargs.get("normval", None)
        titstr = self.config + " " + self.specs

        if qs is None:
            qs = ["1-ME", "2-ME", "3-ME"]
        hstnames = {"1-KE": r"$<\rho v_r^2>$", "2-KE": r"$<\rho v_\theta^2>$", "3-KE": r"$<\rho v_\phi^2>$", "1-ME": r"$<B_r^2>$", "2-ME": r"$<B_\theta^2>$", "3-ME": r"$<B_\phi^2>$", "1-mom":r"$\rho v_r$", "2-mom":r"$\rho v_\theta$", "3-mom":r"$\rho v_\phi$", "tot-E":r"$E$", "mass":r"$M$"}

        # convert time to orbits
        time = time/(2*np.pi)

        plt.figure(); n=0
        if norm:
            if normval is None:
                # find greatest starting value
                starting = [qs[q][0] for q in qs.keys()]
                normval = np.max(np.abs(starting))
            titstr += "\n Normalized to value of {}".format(self.toscinot(normval))

        for q, data in qs.items():
            # normalize data to initial value
            if norm:
                data = data/normval

            plt.plot(time, data, ls=ls[n], label=hstnames[leg_lab[n]])
            n += 1

        plt.xlabel(xlabel)
        if ylabel is not None:
            plt.ylabel(ylabel)
        else:
            plt.ylabel("$<B^2(t)>/<B^2(0)>$")

        if ylog:
            plt.yscale('log')
        plt.legend(loc='best')
        plt.title(titstr)

    #------------------------------------------------------------------------
    def plot_slice(self, q, t, slicedata, grid, slicestr="x2s", H_units=False, sliceind=None, **kwargs):
        """
        t is tind, not tval. tval = self.times[tind, 1]
        options:
        labels
        leg_lab
        vlims: colormap lims
        coord_lims: coordinate limits
        title

        TO DO: XXX
        mix coord_lims/coord_inds: input either coord or index
        time units?
        log

        Can normalize across all time steps via vlim input
        """
        fracnorm = kwargs.get("fracnorm", False)
        levels = kwargs.get("levels", None)
        contourf = kwargs.get("contourf", False)
        diff = kwargs.get("diff", False)
        aAvg = kwargs.get("aAvg", False)

        if len(q.split("-")) > 1 and q.split("-")[1] == "diff":
            diff = True
        f = plt.figure("slice")
        plt.clf(); ax = plt.gca()
        if H_units:
            conversion = 1/self.H
            labend = "$[H]$"
        else:
            conversion = 1.
            labend = "[code-units]"
        # get data
        leg_lab = ""
        [xgrid, ygrid] = grid
        # rescale
        xgrid = xgrid*conversion; ygrid = ygrid*conversion

        if self.logr:
            rstr = r"$\log_{10}(r)$"
        else:
            rstr = r"$r\ $"

        # Fractional
        # make sure no zeros
        # if diff:
            # leg_lab += r"$[$"+self.varnames[q]+r"$(t)-$"+self.varnames[q]+r"$(t=0)]$"
            # if fracnorm:
                # print("fracnorm is true.")
                # if (t0data != 0).all():
                    # slicedata = slicedata/t0data # normalize to initial data
                    # leg_lab += "/"+self.varnames[q]+r"$(t=0)$"
                # else:
                    # print("error: some values in initial data are zero.")
                    # print("not normalizing to initial data.")
        if 'leg_lab' not in kwargs:
            if aAvg:
                leg_lab = "$\langle " + self.varnames[q].strip('$') + r"\rangle_\phi$"
            else:
                leg_lab = self.varnames[q]
        else:
            leg_lab = kwargs['leg_lab']

        if slicestr=="x2s":
            if sliceind is None:
                if self.nx2 % 2 == 0:
                    pos = np.mean(self.x2v[int(self.nx2/2-1):int(self.nx2/2+1)])*conversion
                else:
                    pos = self.xvf[int(self.nx2/2)]*conversion
            else:
                pos = self.x2v[sliceind]*conversion
            if 'labels' not in kwargs:
                labels = [rstr+r"~$\cos\phi$~"+labend, rstr+r"~$\sin\phi$~"+labend]
        elif slicestr=="x3s":
            if sliceind is None:
                if self.nx3 % 2 == 0:
                    pos = np.mean(self.x3v[int(self.nx3/2-1):int(self.nx3/2+1)])*conversion
                else:
                    pos = self.x3v[int(self.nx3/2)]*conversion
            else:
                pos = self.x3v[sliceind]*conversion
            if 'labels' not in kwargs:
                labels = [rstr+r"~$\sin\theta$~"+labend, rstr+r"~$\cos\theta$~"+labend]

        # set defaults
        codef = np.array([-self.x1v[-1], self.x1v[-1], -self.x1v[-1], self.x1v[-1]])*conversion
        coord_lims = kwargs.get("coord_lims", codef)
        if 'title' not in kwargs:
            if slicestr=="x2s":
                title = "Midplane slice"
            elif slicestr=="x3s":
                title = "Vertical slice"

            # title += self.config + "_" + self.specs
            if self.tvals.size <= 1:
                title += ", $t={:05.2f}$ (${:05.2f}$ orbits)".format(self.tvals[0], self.tvals[0]/(2*np.pi))
            else:
                tloc = np.where(self.tinds == t)[0][0]
                title += ", $t={:05.2f}$ (${:05.2f}$ orbits)".format(self.tvals[tloc], self.tvals[tloc]/(2*np.pi))
        else:
            title = kwargs['title']


        if 'vlims' not in kwargs:
            # get extrema
            vlims = [np.nanmin(slicedata), np.nanmax(slicedata)]
            if vlims[0] <= 0:
                # symmetrize
                xcol = "black"
                vlims[1] = np.nanmax([np.abs(vlims[0]), vlims[1]])
                vlims[0] = -vlims[1]
                if "cmap" not in kwargs:
                    cmap = plt.cm.RdBu
                else:
                    cmap = kwargs["cmap"]
                if levels is None:
                    norm = colors.SymLogNorm(linthresh=(vlims[1]-vlims[0])/1000., vmin=vlims[0], vmax=vlims[1])
                else:
                    norm = colors.BoundaryNorm(levels, ncolors=cmap.N)
            else:
                cmap = kwargs.get("cmap", plt.cm.viridis)
                if levels is None:
                    norm = colors.LogNorm(vmin=vlims[0], vmax=vlims[1])
                else:
                    norm = colors.BoundaryNorm(levels, ncolors=cmap.N)
        else:
            vlims = kwargs['vlims']
            if "cmap" not in kwargs:
                if vlims[0] <= 0:
                    cmap = plt.cm.RdBu
                    if levels is None:
                        norm = colors.SymLogNorm(linthresh=(vlims[1]-vlims[0])/1000., vmin=vlims[0], vmax=vlims[1])
                    else:
                        norm = colors.BoundaryNorm(levels, ncolors=cmap.N)
                else:
                    cmap=plt.cm.viridis
                    if levels is None:
                        norm = colors.LogNorm(vmin=vlims[0], vmax=vlims[1])
                    else:
                        norm = colors.BoundaryNorm(levels, ncolors=cmap.N)
            else:
                cmap = kwargs["cmap"]
                if levels is None:
                    norm = colors.Normalize(vmin=vlims[0], vmax=vlims[1])
                else:
                    norm = colors.BoundaryNorm(levels, ncolors=cmap.N)

        if vlims[1] - vlims[0] == 0 or np.isnan(vlims).any():
            print("no data!")
            print("exiting...")
            return 0;

        print("vlims: [{}, {}]".format(vlims[0], vlims[1]))
        # plot!
        if contourf:
            im = ax.contourf(xgrid, ygrid, slicedata, cmap=cmap, levels=levels)
        else:
            im = ax.pcolormesh(xgrid, ygrid, slicedata, cmap=cmap, norm=norm)
        ax.set_xlabel(labels[0])
        ax.set_ylabel(labels[1])
        ax.set_title(title, y=1.1)
        ax.set_xlim((coord_lims[0], coord_lims[1]))
        ax.set_ylim((coord_lims[2], coord_lims[3]))
        ax.set_aspect('equal')
        # ax.axvline(0, ls='--', c='k')
        # ax.axhline(0, ls='--', c='k')
        start, end = ax.get_xlim()
        ax.xaxis.set_ticks(np.linspace(start, end, 5))
        start, end = plt.gca().get_ylim()
        ax.yaxis.set_ticks(np.linspace(start, end, 5))
        plt.xticks(rotation=30)

        cbar = plt.colorbar(im, label=leg_lab, pad=0.1)
        plt.tight_layout()


    #------------------------------------------------------------------------
    def plot_2slice(self, q, t, data, grids, H_units=False, **kwargs):
        """
        options:
        labels
        leg_lab
        vlims: colormap lims
        coord_lims: coordinate limits
        title

        TO DO: XXX
        mix coord_lims/coord_inds: input either coord or index
        time units?
        log

        Can normalize across all time steps via vlim input
        """
        [data_mp, data_vert] = data
        fracnorm = kwargs.get("fracnorm", False)
        aAvg = kwargs.get("aAvg", False)
        leg_lab = ""
        if H_units:
            conversion = 1/self.H
        else:
            conversion = 1.
        # set defaults
        diff = kwargs.get("diff", False)
        if 'title' not in kwargs:
            if self.times is not None:
                if self.tvals.size <= 1:
                    title = "$t={:05.2f}$ (${:05.2f}$ orbits)".format(self.tvals[0], self.tvals[0]/(2*np.pi))
                else:
                    tloc = np.where(self.tinds == t)[0][0]
                    title = "$t={:05.2f}$ (${:05.2f}$ orbits)".format(self.tvals[tloc], self.tvals[tloc]/(2*np.pi))
            else:
                title += "$t={:d}$".format(t)
        else:
            title = kwargs['title']
        codef = np.array([-self.x1v[-1], self.x1v[-1], -self.x1v[-1], self.x1v[-1]])*conversion
        coord_lims = kwargs.get("coord_lims", codef)
        if self.logr:
            rstr = r"$\log_{10}(r)$"
        else:
            rstr = r"$r\ $"

        # get data
        [x2s_grid, x3s_grid] = grids
        # Fractional
        # make sure no zeros
        # if diff:
            # if 'leg_lab' not in kwargs:
                # leg_lab += r"$[$"+self.varnames[q]+r"$(t)-$"+self.varnames[q]+r"$(t=0)]$"
            # if fracnorm:
                # if (t0data != 0).all():
                    # raw_data = raw_data/t0data # normalize to initial data
                    # leg_lab += "/" + self.varnames[q]+r"$(t=0)$"
                # else:
                    # print("error: some values in initial data are zero.")
        if 'leg_lab' not in kwargs:
            if aAvg:
                leg_lab = "$\langle " + self.varnames[q].strip('$') + r"\rangle_\phi$"
            else:
                leg_lab = self.varnames[q]
        else:
            leg_lab = kwargs['leg_lab']

        nr = 1; nc = 2; fign = 1
        if 'vlims' not in kwargs:
            # get extrema
            vmin = np.nanmin([np.nanmin(data_vert), np.nanmin(data_mp)])
            vmax = np.nanmax([np.nanmax(data_vert), np.nanmax(data_mp), np.nanmax(data_mp)])
            vlims = [vmin, vmax]
        else:
            vlims = kwargs['vlims']
        if vlims[1] - vlims[0] == 0 or np.isnan(vlims).any():
            print("no data.")
            print("Exiting...")
            return 0;

        if 'vlims' not in kwargs and vlims[0] <= 0:
            # Symmetrize
            vlims[1] = np.max([np.abs(vlims[0]), vlims[1]])
            vlims[0] = -vlims[1]
            norm = colors.SymLogNorm(linthresh=(vlims[1]-vlims[0])/1000., vmin=vlims[0], vmax=vlims[1])
            xcol = "black"
            if "cmap" not in kwargs:
                cmap = plt.cm.RdBu
            else:
                cmap = kwargs["cmap"]
        else:
            norm = colors.LogNorm(vmin=vlims[0], vmax=vlims[1])
            if "cmap" not in kwargs:
                if vlims[0] < 0:
                    cmap = plt.cm.RdBu
                else:
                    cmap = plt.cm.viridis
            else:
                cmap = kwargs["cmap"]

        fign = plt.figure(1, figsize=(30,50))
        plt.clf()
        print("vlims: [{}, {}]".format(vlims[0], vlims[1]))

        grid = AxesGrid(fign, 111, nrows_ncols=(nr, nc), cbar_location="right", cbar_pad=0.05,
                        cbar_mode="single", label_mode="all", axes_pad=1.0, aspect=True)

        imy = grid[0].pcolormesh(x2s_grid[0]*conversion, x2s_grid[1]*conversion, data_mp, cmap=cmap, norm=norm, vmin=vlims[0], vmax=vlims[1])
        imz = grid[1].pcolormesh(x3s_grid[0]*conversion, x3s_grid[1]*conversion, data_vert, cmap=cmap, norm=norm, vmin=vlims[0], vmax=vlims[1])

        plt.subplots_adjust(top=0.80) # get suptitle in the right spot
        grid[1].set_title("Vertical slice")
        grid[0].set_title("Midplane slice")
        if H_units:
            labend = "$[H]$"
        else:
            labend = "[code-units]"

        grid[0].set_xlabel(rstr+r"~$\cos\phi$~"+labend)
        grid[0].set_ylabel(rstr+r"~$\sin\phi$~"+labend)
        grid[1].set_xlabel(rstr+r"~$\cos\theta$~"+labend)
        grid[1].set_ylabel(rstr+r"~$\sin\theta$~"+labend)
        grid[0].set_xlim((coord_lims[0], coord_lims[1]))
        grid[0].set_ylim((coord_lims[0], coord_lims[1]))
        grid[1].set_xlim((coord_lims[2], coord_lims[3]))
        grid[1].set_ylim((coord_lims[2], coord_lims[3]))
        # grid[0].axvline(0, ls='--', c='k')
        # grid[0].axhline(0, ls='--', c='k')
        grid[0].tick_params(axis='x', rotation=30)
        # grid[1].axvline(0, ls='--', c='k')
        # grid[1].axhline(0, ls='--', c='k')
        grid[1].tick_params(axis='x', rotation=30)

        # ticks = np.linspace(vlims[0], vlims[1], 7)
        # ticklabs = [self.toscinot(val) for val in ticks]
        # grid.cbar_axes[0].colorbar(imy, format=LogFormatterSciNotation(), ticks=ticks);
        grid.cbar_axes[0].colorbar(imy) #, format=LogFormatterSciNotation());
        # grid.cbar_axes[0].set_yticklabels(ticklabs)
        grid.cbar_axes[0].set_ylabel(leg_lab)

        plt.subplots_adjust(top=1.00) # get suptitle in the right spot
        plt.suptitle(title, horizontalalignment='center', verticalalignment='top', fontsize=15)

        return fign

    #------------------------------------------------------------------------

    def plot_profile(self, q, t, slicedata, sstr, sind=None, **kwargs):
        logy = kwargs.get("logy", False)
        logx = kwargs.get("logx", False)
        title = kwargs.get("title", None)
        xlab = kwargs.get("xlab", None)
        # XX xlim different for radial vs vertical profile
        # xlim = kwargs.get("xlim", [self.x1v[0], self.x1v[-1]])
        ylim = kwargs.get("ylim", None)
        toavg = kwargs.get("avg", False) # note: is only for x2s midplane
        if title is None:
            title = self.config + "-" + self.specs + "\n"
        ylab = kwargs.get("ylab", None)
        if ylab is None:
            if toavg:
                ylab = "$\langle " + self.varnames[q].strip('$') + r"\rangle_\phi$"
            else:
                ylab = self.varnames[q]

        if sstr == "radial": # midplane radial profile
            # sind is phi position...can only do if have vertical slice there and for non-azimuthally-averaged
            if sind is None and not toavg:
                sind = 0
            conversion = kwargs.get("conversion", 1.) 
            xvals = self.x1v
            # average two surround cells if necessary
            if self.nx2 % 2 == 0:
                ydata = np.mean([slicedata[int(self.nx2/2-1), :int(self.nx1)], slicedata[int(self.nx2/2+1), :int(self.nx1)]], axis=0)
            else:
                ydata = slicedata[int(self.nx2/2), :int(self.nx3)]

            xlim = kwargs.get("xlim", [xvals[0], xvals[-1]])
            if xlim is None:
                xlim = [xvals[0], xvals[-1]]
            xlim.sort()
            mask = (xvals >= xlim[0]) & (xvals <= xlim[-1])
            xvals = xvals[mask]
            ydata = ydata[mask]

            if xlab is None:
                xlab = "$r/R_0$"
            if title not in kwargs:
                if toavg:
                    title += "Radial profile averaged over $\phi$"
                else:
                    title += "Radial profile at $\phi={:05.2f}$".format(self.x3v[sind])
        elif sstr == "vertical": # vertical profile
            # sind is r position
            if sind is None:
                # find index such that r = 1.0
                sind = (np.abs(data_plot.x1v - 1.0)).argmin()
            conversion = kwargs.get("conversion", 1/self.x1v[sind]) # H(r=1)
            xvals = self.x3s_grid[1][:int(self.nx2), sind]*conversion
            xlim = kwargs.get("xlim", [xvals[0], xvals[-1]])
            if xlim is None:
                xlim = [xvals[0], xvals[-1]]
            xlim.sort()

            # slice data
            ydata = slicedata[:int(self.nx2), sind]
            # only use data within xlim:
            mask = (xvals >= xlim[0]) & (xvals <= xlim[-1])
            xvals = xvals[mask]
            ydata = ydata[mask]

            if xlab is None:
                xlab = r"$z/R$"
            if "title" not in kwargs:
                title += r"Vertical (along $\theta$) profile at $r={:05.2f}$".format(self.x1v[sind])
        else:
            print("Error: unknown string option " + sstr)
            print("Available options: x2s or x3s.")
            return

        if "title" not in kwargs:
            tloc = np.where(self.tinds == t)[0][0]
            title += "\n Time: ${:05.2f}$ (${:05.2f}$ orbits)".format(self.tvals[tloc], self.tvals[tloc]/(2*np.pi))

        plt.figure()
        plt.plot(xvals, ydata, marker='o', markersize=2) # XX linestyle
        if logy:
            plt.yscale("log")
        else:
            plt.gca().axhline(0, color='black', ls='--')
        if logx:
            plt.xscale("log")
        if sstr == "x2s":
            try:
                plt.gca().axvline(self.r1, color='black', ls='--')
                plt.gca().axvline(self.r2, color='black', ls='--')
                plt.gca().axvline(self.rm, color='black', ls='--')
            except: pass
        plt.title(title)
        plt.xlabel(xlab)
        plt.ylabel(ylab)
        plt.gca().set_xlim(xlim)
        if ylim is not None:
            plt.gca().set_ylim(ylim)
        plt.tight_layout()

    #------------------------------------------------------------------------

    def fig_save(self, **kwargs):
        import datetime
        from PIL import Image, PngImagePlugin
        metadata={"Creator":"Lia Hankla", "Description":"Simulation: "+self.config+"_"+self.specs,
            "Creation Time":datetime.datetime.now().isoformat()}
        if "fig_type" not in kwargs:
            fig_type = "scratch/"
        else:
            fig_type = kwargs["fig_type"]
        if "fig_name" not in kwargs:
            tit = plt.gca().get_title()
            tit = tit.replace(' ', '-')
            fig_name = tit + ".png"
        else:
            fig_name = kwargs["fig_name"]
        fig_dir = self.fig_save_fp + fig_type + "/"

        if not os.path.exists(fig_dir):
            os.makedirs(fig_dir)
        fig_name = fig_dir + fig_name
        try:
            plt.savefig(fig_name, bbox_inches='tight')
            # insert metadata after the fact
            img = Image.open(fig_name)
            meta = PngImagePlugin.PngInfo()
            for k, v in metadata.items():
                meta.add_text(k, v)
            img.save(fig_name, "png", pnginfo=meta)
        except:
            print("Could not save figure " + fig_name)
        else:
            print("saving figure " + fig_name)



    # -----------------------------
    def toscinot(self, val):
        if val == 0.:
            return r"${:.03f}$".format(0.)
        exp = int(np.log10(np.abs(val)))
        if exp < 0:
            exp = exp - 1
        base = val/10.**exp
        if base == 10:
            base = 1.; exp = exp + 1
        if exp == 0:
            return r"${:.03f}$".format(base, exp)
        else:
            return r"${:.03f}\times10^{{{:d}}}$".format(base, exp)
