# Mac_Low_Turbulence_Driving
A C++ code using [Mac Low(1999)](https://ui.adsabs.harvard.edu/abs/1999ApJ...524..169M/abstract)'s method to drive turbulence, with an ORION2-format output.

## How to compile and use the code:
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

After running the program, just follow the program's prompts to set the number of grids, wavenumber, and random seeds. Please note that the wavenumber you set should be the maximum wavenumber. For example, if you expect the power spectrum to be a top-hat function of \[1,2\], you should set the wave number to 2.

# Warning: Still testing! Please use with caution and verify results.
