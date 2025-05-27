# Mac_Low_Turbulence_Driving
A C++ code using [Mac Low(1999)](https://ui.adsabs.harvard.edu/abs/1999ApJ...524..169M/abstract)'s method to drive turbulence, with an ORION2-format output.

## Compile and use the code:
**Dependencies**: This code relies on FFTW3 and HDF5 libraries. You need to install them in your system first. For Mac Silicon users, homebrew is highly recommended, the library path is `/opt/homebrew/lib` while the include path is `/opt/homebrew/include`. For Linux users, my personal advice is that you can use `--prefix` to specify the install path for every libraries (not only for using this code), because not every code can find the dependency libraries (especially hdf5) installed in the default location well. Moreover, for numerical simulation workers, sometimes you need to use different configurations of dependency libraries, and you may keep several different configurations of the same dependency library at the same time. At this time, using the default location will cause confusion.

**Installation**:

**Step 1**: Get the code and enter the directory. 
```bash
git clone https://github.com/tianninglyu/Mac_Low_Turbulence_Driving.git`
cd Mac_Low_Turbulence_Driving
```

**Step 2**: Compile the code. 
```bash
g++ -o TurbDrive main.cpp TurbDrive.cpp -I/path/to/include -L/path/to/lib
```
Now you can see a executable file called 'TurbDrive' in the working directory.

**Step 3**: Run the code.

After running the program, just follow the program's prompts to set the number of grids, wavenumber, and random seeds. Please note that the wavenumber you set should be the maximum wavenumber. For example, if you expect the power spectrum to be a top-hat function of \[1,2\], you should set the wavenumber to 2.

You will finally get a perturbation field file named 'zdrv_n\[number\]_k\[number\].hdf5'. For example, if `n=128` and `k=2`, the filename will be `zdrv_n128_k2.hdf5`. For using it in ORION2, you should rename it to `zdrv.hdf5`.

## Citation
This code is just an implementation of the turbulence driven method proposed by [Mac Low(1999)](https://ui.adsabs.harvard.edu/abs/1999ApJ...524..169M/abstract). If you use this code in academic publications, you should cite [Mac Low(1999)](https://ui.adsabs.harvard.edu/abs/1999ApJ...524..169M/abstract). 
```latex
@ARTICLE{1999ApJ...524..169M,
       author = {{Mac Low}, Mordecai-Mark},
        title = "{The Energy Dissipation Rate of Supersonic, Magnetohydrodynamic Turbulence in Molecular Clouds}",
      journal = {\apj},
     keywords = {ISM: CLOUDS, ISM: KINEMATICS AND DYNAMICS, ISM: MAGNETIC FIELDS, MAGNETOHYDRODYNAMICS: MHD, TURBULENCE, ISM: Clouds, ISM: Kinematics and Dynamics, ISM: Magnetic Fields, Magnetohydrodynamics: MHD, Turbulence, Astrophysics},
         year = 1999,
        month = oct,
       volume = {524},
       number = {1},
        pages = {169-178},
          doi = {10.1086/307784},
archivePrefix = {arXiv},
       eprint = {astro-ph/9809177},
 primaryClass = {astro-ph},
       adsurl = {https://ui.adsabs.harvard.edu/abs/1999ApJ...524..169M},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```
I also kindly request you to mention this code and include the GitHub URL in the footnotes or acknowledgments section of your publication.

# Warning: Still testing! Please use with caution and verify results.
