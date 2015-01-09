# This is a wrapper for JAMM to use it in Galaxy
# It will be accompanied with a jamm.xml file, which specifies the interface that Galaxy is showing to the user
# as well as the JAMM software including JAMM.sh, the bash script that is actually called.

# This wrapper does the following things:
# map the files (from Galaxy history) to a directory
# pass the parameters from the GUI to JAMM.sh
# call JAMM.sh
# maybe map the resulting tabular files back to history items, the XML could take care of that if output is fixed

#import optparse
import argparse, os, shutil, subprocess, sys, tempfile
import shlex
# importing some of the modules used, especially for symlinking, dir manipulation and tool calling.
# since python2.7, argparse is preferred over optparse.

# for reference wrappers, see
# https://bitbucket.org/fubar/rossdev/src/4a91f99b5e1953270c9b35d2ca70c325a10fcfcb/rgdev/bwa_wrapper.py?at=default
# https://www.e-biogenouest.org/wiki/WrapperpythonforGalaxy:Asyntaxexample

def main(): 

    #Read command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', dest = 'input', nargs='+')
    parser.add_argument('-g', dest = 'gsize')
    parser.add_argument('-z', dest = 'peakfile')
    
    parser.add_argument('-m', dest = 'mode')
    parser.add_argument('-r', dest = 'resolution')
    parser.add_argument('-p', dest = 'processes')
    #parser.add_argument('-b', dest = 'binSize')
    
    args = parser.parse_args()
    
    print "################################################" 
    print "Wrapper debugging" 
    print "################################################" 
    print "Files to be used:"
    for j in args.input:
        print j

    print "output file:"
    print args.peakfile
    print "current working dir:"
    print os.getcwd() 
    print "dir with jammwrapper in it:"
    path = os.path.abspath(os.path.dirname(sys.argv[0]))
    print path

    # optparse was depracted, can still be found in may of the example wrappers
    #parser = optparse.OptionParser()
    #parser.add_option( '-i',  dest='input', help='input bed files' )
    #parser.add_option( '-c',  dest='csize', help='chr sizes' )
    #(options, args) = parser.parse_args()	
   
    # create temp dir
    tmp_dir = tempfile.mkdtemp()
    #os.chdir(tmp_dir)
    # symlink creation
    for file in args.input:
        filen =  tmp_dir + "/" + os.path.basename(os.path.splitext(file)[0])+".bed"
        os.symlink(file, filen)
        
        # in case temp files should have random names
        # ref_file = tempfile.NamedTemporaryFile( dir=tmp_dir )
        # ref_file_name = ref_file.name
        # ref_file.close()
        # os.symlink( file, ref_file_name )

	    
    #command = ( "/home/cmesser/jamm/JAMM.sh -s %s -g %s -o results -m %s -r %s -p %s"
    command = ( "bash %s/JAMM.sh -s %s -g %s -o results -m %s -r %s -p %s"
     % ( path, tmp_dir, args.gsize, args.mode, args.resolution, args.processes ) ) 
    print command
    # depending on how your programm is called, it may be necessary to use shlex.split on the command string before
    # in this case, this was actually harmful. idk why
#    command = shlex.split(command)
    # Please read the docs for subprocess.Popen. You can't just pass as string as variables need to be extended
    # Note that shell = True can be a security hazard. 
#    proc = subprocess.Popen( command, executable = '/bin/bash', shell=True )
    proc = subprocess.Popen( command, shell=True)
    returncode = proc.wait()
    
    #mv files to a place where galaxy wrapper can find them
    # mvcommand = "mv %s/results/peaks/all.peaks.narrowPeak %s" % ( tmp_dir, args.peakfile )
    mvcommand = "mv results/peaks/all.peaks.narrowPeak %s" % args.peakfile
    
    print mvcommand
    os.system(mvcommand)

# clean up temp dir
    if os.path.exists( tmp_dir ):
        shutil.rmtree( tmp_dir )    

    
if __name__ == "__main__":
    main()

