#include "langevinDynamicsDoubleBarrier.h"

using namespace std;
using namespace LANGEVINEQUATION;

int main(int argc, char* argv[]){

    MPI_Init(&argc, &argv);
    int nodeIndex = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &nodeIndex);

    // ArgumentParser inputParser;
    // inputParser.AddOptionalArgument("-o", "--output", "Output", 1, {"dynamicBarrier"});
    // inputParser.AddOptionalArgument("-x0", "--x0", "X0", 1, {"0.0"});
    // inputParser.AddOptionalArgument("-v0", "--v0", "V0", 1, {"1.0"});
    // inputParser.AddOptionalArgument("-g", "--gamma", "Gamma", 1, {"1.0"});
    // inputParser.AddOptionalArgument("-kb", "--kb", "Boltzman factor", 1, {"1.0"});
    // inputParser.AddOptionalArgument("-t", "--temperature", "Temperature", 1, {"1.0"});
    // inputParser.AddOptionalArgument("-m", "--mass", "Mass", 1, {"1.0"});
    // inputParser.AddOptionalArgument("-d", "--dt", "Timestep", 1, {"0.1"});
    // inputParser.AddOptionalArgument("-l", "--length", "Number of steps", 1, {"10000"});
    // inputParser.AddOptionalArgument("-f1", "--savefreq", "Save frequency", 1, {"1"});
    // inputParser.AddOptionalArgument("-f0", "--statfreq", "Statistics frequency", 1, {"1"});
    // inputParser.AddOptionalArgument("-vr", "--vrandom", "Randomize initial velocity", 1, {"true"});
    // inputParser.AddOptionalArgument("-er", "--etarandom", "Randomize initial eta", 1, {"true"});
    // inputParser.AddOptionalArgument("-sr", "--staterandom", "Randomize initial state", 1, {"true"});
    // inputParser.AddOptionalArgument("-ss", "--savestate", "Same state of barriers", 1, {"false"});
    // inputParser.AddOptionalArgument("-p", "--periodic", "Periodic length", 1, {"10.0"});
    // inputParser.AddOptionalArgument("-b0", "--barrier0", "Barrier height 0", 2, {"2.0", "6.0"});
    // inputParser.AddOptionalArgument("-b1", "--barrier1", "Barrier height 1", 2, {"2.0", "6.0"});
    // inputParser.AddOptionalArgument("-s0", "--initialstate", "Initial state", 1, {"0"});
    // inputParser.AddOptionalArgument("-sl0", "--slopehalflength0", "Slope half length 0", 1, {"1.0"});
    // inputParser.AddOptionalArgument("-sl1", "--slopehalflength1", "Slope half length 1", 1, {"1.0"});
    // inputParser.AddOptionalArgument("-bp0", "--barrierposition0", "Barrier position 0", 1, {"1.0"});
    // inputParser.AddOptionalArgument("-bp1", "--barrierposition1", "Barrier position 1", 1, {"2.0"});
    // inputParser.AddOptionalArgument("-q", "--q", "Q", 1, {"0.2"});
    // inputParser.AddOptionalArgument("-tau", "--tau", "Tau", 1, {"0.1"});
    // inputParser.AddOptionalArgument("-eta0", "--eta0", "Eta0", 1, {"0.0"});
    // inputParser.AddOptionalArgument("-wt", "--waitingtime", "Waiting time", 1, {"false"});
    // inputParser.AddOptionalArgument("-fpt", "--fpt", "First passage time", 1, {"false"});
    // inputParser.AddOptionalArgument("-tt", "--tpt", "transition path time", 1, {"false"});
    // inputParser.AddOptionalArgument("-fptcg", "--fptcg", "First passage time CG", 2, {"false", "10"}, false, ArgumentAction::Push);
    // inputParser.AddOptionalArgument("-rr", "--readrestart", "Read restart", 1, {"false"});
    // inputParser.AddOptionalArgument("-wr", "--writerestart", "Write restart", 1, {"false"});
    // inputParser.AddOptionalArgument("-ir", "--inputrestart", "Input restart file base", 1, {"inputRestart"});    
    // inputParser.AddOptionalArgument("-or", "--outputrestart", "Output restart file base", 1, {"outputRestart"});
    // inputParser.AddOptionalArgument("-r", "--run", "Number of runs", 1, {"100000"});
    // inputParser.AddOptionalArgument("-pl", "--previouslength", "Number of previous steps", 1, {"0"});
    // inputParser.AddOptionalArgument("-we", "--writeeta", "Write eta", 1, {"false"});
    // inputParser.AddOptionalArgument("-pm", "--probmin", "Prob hist min", 1, {"-100.0"});
    // inputParser.AddOptionalArgument("-pa", "--probmax", "Prob hist max", 1, {"100.0"});
    // inputParser.AddOptionalArgument("-pb", "--probbin", "Prob hist bin", 1, {"0.1"});
    // inputParser.AddOptionalArgument("-pt", "--probtime", "Time of prob hist", 1, {"0.1"}, true, ArgumentAction::Push);
    // inputParser.AddOptionalArgument("-wx", "--writex", "Write x", 1, {"false"});
    // inputParser.AddOptionalArgument("-sp", "--sample", "Sample per order", 1, {"20"});
    // inputParser.AddOptionalArgument("-rng", "--random", "Random number seed", 1, {"-1"});
    // inputParser.AddOptionalArgument("-xcr", "--xcorrelation", "Coordinate correlation", 1, {"false"});
    // inputParser.AddOptionalArgument("-bt", "--btimer", "Print timer", 1, {"false"});
    // inputParser.AddOptionalArgument("-v0d", "--v0d", "V0 distribution", 1, {"false"});
    // inputParser.AddOptionalArgument("-eta0d", "--eta0d", "Eta0 distribution", 1, {"false"});
    // inputParser.AddOptionalArgument("-s0d", "--s0d", "State0 distribution", 1, {"false"});
    // inputParser.AddOptionalArgument("-od", "--overdamped", "Overdamped", 1, {"false"});
    // inputParser.AddOptionalArgument("-fh", "--fpthist", "FPT hist min max nbin", 3, {"0.0", "10.0", "100"});
    // inputParser.AddOptionalArgument("-eq", "--eq", "Equal prob, K12", 2, {"true", "1.0"});

    // inputParser.Parse(&argc, &argv);

    INPUTPARSER::InputParser inputParser;
    inputParser.AddOption(/* -o    */ "OUTPUTBASE", "Output base", "langevinDynamicsSingleBarrier", INPUTPARSER::InputNumberRequirement::Equal, 1, INPUTPARSER::InputOptionType::MustAndOnce, INPUTPARSER::InputType::Output);
    inputParser.AddOption(/* -od   */ "OVERDAMPED", "Whether use overdamped equation", "false", INPUTPARSER::InputNumberRequirement::Equal, 1, INPUTPARSER::InputOptionType::MustAndOnce, INPUTPARSER::InputType::Parameter);
    inputParser.AddOption(/* -x0   */ "X0", "Initial position (x0)", "0.0", INPUTPARSER::InputNumberRequirement::Equal, 1, INPUTPARSER::InputOptionType::MustAndOnce, INPUTPARSER::InputType::Parameter);
    inputParser.AddOption(/* -v0   */ "V0", "Initial velocity (v0)", "1.0", INPUTPARSER::InputNumberRequirement::Equal, 1, INPUTPARSER::InputOptionType::MustAndOnce, INPUTPARSER::InputType::Parameter);
    inputParser.AddOption(/* -g    */ "GAMMA", "Friction coefficient (Gamma)", "1.0", INPUTPARSER::InputNumberRequirement::Equal, 1, INPUTPARSER::InputOptionType::MustAndOnce, INPUTPARSER::InputType::Parameter);
    inputParser.AddOption(/* -kb   */ "KB", "Boltzman factor", "1.0", INPUTPARSER::InputNumberRequirement::Equal, 1, INPUTPARSER::InputOptionType::MustAndOnce, INPUTPARSER::InputType::Parameter);
    inputParser.AddOption(/* -t    */ "TEMPERATURE", "Temperature", "1.0", INPUTPARSER::InputNumberRequirement::Equal, 1, INPUTPARSER::InputOptionType::MustAndOnce, INPUTPARSER::InputType::Parameter);
    inputParser.AddOption(/* -m    */ "MASS", "Mass", "1.0", INPUTPARSER::InputNumberRequirement::Equal, 1, INPUTPARSER::InputOptionType::MustAndOnce, INPUTPARSER::InputType::Parameter);
    inputParser.AddOption(/* -vr   */ "RANDOMIZEV", "Whether randomize initial velocity", "true", INPUTPARSER::InputNumberRequirement::Equal, 1, INPUTPARSER::InputOptionType::MustAndOnce, INPUTPARSER::InputType::Parameter);

    inputParser.AddOption(/* -p    */ "PERIODICLENGTH", "Periodic length", "10.0", INPUTPARSER::InputNumberRequirement::Equal, 1, INPUTPARSER::InputOptionType::MustAndOnce, INPUTPARSER::InputType::Parameter);
    inputParser.AddOption(/* -b0   */ "BARRIERHEIGHTS0", "Barrier heights 0 (two)", "2.0 6.0", INPUTPARSER::InputNumberRequirement::Equal, 2, INPUTPARSER::InputOptionType::MustAndOnce, INPUTPARSER::InputType::Parameter);
    inputParser.AddOption(/* -b1   */ "BARRIERHEIGHTS1", "Barrier heights 1 (two)", "2.0 6.0", INPUTPARSER::InputNumberRequirement::Equal, 2, INPUTPARSER::InputOptionType::MustAndOnce, INPUTPARSER::InputType::Parameter);
    inputParser.AddOption(/* -sl0  */ "BARRIERHALFWIDTH0", "Barrier half width 0", "1.0", INPUTPARSER::InputNumberRequirement::Equal, 1, INPUTPARSER::InputOptionType::MustAndOnce, INPUTPARSER::InputType::Parameter);
    inputParser.AddOption(/* -sl1  */ "BARRIERHALFWIDTH1", "Barrier half width 1", "1.0", INPUTPARSER::InputNumberRequirement::Equal, 1, INPUTPARSER::InputOptionType::MustAndOnce, INPUTPARSER::InputType::Parameter);
    inputParser.AddOption(/* -bp0  */ "BARRIERPOSITION0", "Barrier position 0", "1.0", INPUTPARSER::InputNumberRequirement::Equal, 1, INPUTPARSER::InputOptionType::MustAndOnce, INPUTPARSER::InputType::Parameter);
    inputParser.AddOption(/* -bp1  */ "BARRIERPOSITION1", "Barrier position 1", "2.0", INPUTPARSER::InputNumberRequirement::Equal, 1, INPUTPARSER::InputOptionType::MustAndOnce, INPUTPARSER::InputType::Parameter);

    inputParser.AddOption(/* -d    */ "TIMESTEP", "Timestep (dt)", "0.1", INPUTPARSER::InputNumberRequirement::Equal, 1, INPUTPARSER::InputOptionType::MustAndOnce, INPUTPARSER::InputType::Parameter);
    inputParser.AddOption(/* -l    */ "NUMBEROFSTEPS", "Number of steps", "10000", INPUTPARSER::InputNumberRequirement::Equal, 1, INPUTPARSER::InputOptionType::MustAndOnce, INPUTPARSER::InputType::Parameter);
    inputParser.AddOption(/* -pl   */ "NUMBEROFPREVIOUSSTEPS", "Number of previous steps", "0", INPUTPARSER::InputNumberRequirement::Equal, 1, INPUTPARSER::InputOptionType::MustAndOnce, INPUTPARSER::InputType::Parameter);
    inputParser.AddOption(/* -r    */ "NUMBEROFRUNS", "Number of runs", "1", INPUTPARSER::InputNumberRequirement::Equal, 1, INPUTPARSER::InputOptionType::MustAndOnce, INPUTPARSER::InputType::Parameter);

    inputParser.AddOption(/* -f1   */ "FREQSAVE", "Save frequency", "1", INPUTPARSER::InputNumberRequirement::Equal, 1, INPUTPARSER::InputOptionType::MustAndOnce, INPUTPARSER::InputType::Parameter);
    inputParser.AddOption(/* -f0   */ "FREQSTAT", "Statistics frequency", "1", INPUTPARSER::InputNumberRequirement::Equal, 1, INPUTPARSER::InputOptionType::MustAndOnce, INPUTPARSER::InputType::Parameter);
    inputParser.AddOption(/* -v0d  */ "DISTRIBUTIONV0", "Whether calculate V0 distribution", "false", INPUTPARSER::InputNumberRequirement::Equal, 1, INPUTPARSER::InputOptionType::MustAndOnce, INPUTPARSER::InputType::Parameter);
    inputParser.AddOption(/* -s0d  */ "DISTRIBUTIONS0", "Whether calculate state0 distribution", "false", INPUTPARSER::InputNumberRequirement::Equal, 1, INPUTPARSER::InputOptionType::MustAndOnce, INPUTPARSER::InputType::Parameter);

    inputParser.AddOption(/* -tau  */ "TAU", "Tau of bichotomic noise eta", "0.1", INPUTPARSER::InputNumberRequirement::Equal, 1, INPUTPARSER::InputOptionType::MustAndOnce, INPUTPARSER::InputType::Parameter);
    inputParser.AddOption(/* -s0   */ "INITIALSTATE", "Initial barrier state (0 or 1)", "0", INPUTPARSER::InputNumberRequirement::Equal, 1, INPUTPARSER::InputOptionType::MustAndOnce, INPUTPARSER::InputType::Parameter);
    inputParser.AddOption(/* -sr   */ "RANDOMIZESTATES", "Whether randomize initial states", "true", INPUTPARSER::InputNumberRequirement::Equal, 1, INPUTPARSER::InputOptionType::MustAndOnce, INPUTPARSER::InputType::Parameter);
    inputParser.AddOption(/* -ss   */ "UNIFORMINITALSTATES", "Uniform initial states", "false", INPUTPARSER::InputNumberRequirement::Equal, 1, INPUTPARSER::InputOptionType::MustAndOnce, INPUTPARSER::InputType::Parameter);
    inputParser.AddOption(/* -eq   */ "EQ", "Whether equal prob, constant forward K", "true 1.0", INPUTPARSER::InputNumberRequirement::Equal, 2, INPUTPARSER::InputOptionType::MustAndOnce, INPUTPARSER::InputType::Parameter);

    inputParser.AddOption(/* -wt   */ "WAITINGTIME", "Whether calculate waiting time", "false", INPUTPARSER::InputNumberRequirement::Equal, 1, INPUTPARSER::InputOptionType::MustAndOnce, INPUTPARSER::InputType::Parameter);
    inputParser.AddOption(/* -fpt  */ "FPT", "Whether calculate first passage time", "false", INPUTPARSER::InputNumberRequirement::Equal, 1, INPUTPARSER::InputOptionType::MustAndOnce, INPUTPARSER::InputType::Parameter);
    inputParser.AddOption(/* -tt   */ "TPT", "Whether calculate transition path time", "false", INPUTPARSER::InputNumberRequirement::Equal, 1, INPUTPARSER::InputOptionType::MustAndOnce, INPUTPARSER::InputType::Parameter);
    inputParser.AddOption(/* -fptcg*/ "FPTCG", "Whether calculate coarse grained first passage time", "false 10", INPUTPARSER::InputNumberRequirement::Equal, 2, INPUTPARSER::InputOptionType::Multiple, INPUTPARSER::InputType::Parameter);

    inputParser.AddOption(/* -rr   */ "READRESTART", "Whether read restart", "false", INPUTPARSER::InputNumberRequirement::Equal, 1, INPUTPARSER::InputOptionType::MustAndOnce, INPUTPARSER::InputType::Parameter);
    inputParser.AddOption(/* -wr   */ "WRITERESTART", "Whether write restart", "false", INPUTPARSER::InputNumberRequirement::Equal, 1, INPUTPARSER::InputOptionType::MustAndOnce, INPUTPARSER::InputType::Parameter);
    inputParser.AddOption(/* -ir   */ "INPUTRESTART", "Input restart file base", "inputRestart", INPUTPARSER::InputNumberRequirement::Equal, 1, INPUTPARSER::InputOptionType::MustAndOnce, INPUTPARSER::InputType::Input);
    inputParser.AddOption(/* -or   */ "OUTPUTRESTART", "Output restart file base", "outputRestart", INPUTPARSER::InputNumberRequirement::Equal, 1, INPUTPARSER::InputOptionType::MustAndOnce, INPUTPARSER::InputType::Output);
    inputParser.AddOption(/* -wx   */ "WRITEX", "Whether write x trajectory", "false", INPUTPARSER::InputNumberRequirement::Equal, 1, INPUTPARSER::InputOptionType::MustAndOnce, INPUTPARSER::InputType::Parameter);

    inputParser.AddOption(/* -sp   */ "NUMBEROFLOGSAMPLES", "Number of samples point per order", "20", INPUTPARSER::InputNumberRequirement::Equal, 1, INPUTPARSER::InputOptionType::MustAndOnce, INPUTPARSER::InputType::Parameter);
    inputParser.AddOption(/* -rng  */ "RANDOMSEED", "Random number seed", "-1", INPUTPARSER::InputNumberRequirement::Equal, 1, INPUTPARSER::InputOptionType::MustAndOnce, INPUTPARSER::InputType::Parameter);
    inputParser.AddOption(/* -bt   */ "PRINTTIMER", "Whether print timer", "false", INPUTPARSER::InputNumberRequirement::Equal, 1, INPUTPARSER::InputOptionType::MustAndOnce, INPUTPARSER::InputType::Parameter);
    // inputParser.AddOption(/* -er   */ "RANDOMETA", "Randomize initial eta", "true", INPUTPARSER::InputNumberRequirement::Equal, 1, INPUTPARSER::InputOptionType::MustAndOnce, INPUTPARSER::InputType::Parameter);
    // inputParser.AddOption(/* -q  */ "Q", "Q", "0.2", INPUTPARSER::InputNumberRequirement::Equal, 1, INPUTPARSER::InputOptionType::MustAndOnce, INPUTPARSER::InputType::Parameter);
    // inputParser.AddOption(/* -eta0*/ "ETA0", "Initial state", "0.0", INPUTPARSER::InputNumberRequirement::Equal, 1, INPUTPARSER::InputOptionType::MustAndOnce, INPUTPARSER::InputType::Parameter);
    // inputParser.AddOption(/* -we   */ "WRITEETA", "Write eta", "false", INPUTPARSER::InputNumberRequirement::Equal, 1, INPUTPARSER::InputOptionType::MustAndOnce, INPUTPARSER::InputType::Parameter);
    // inputParser.AddOption(/* -pm   */ "PROBMIN", "Prob hist min", "-100.0", INPUTPARSER::InputNumberRequirement::Equal, 1, INPUTPARSER::InputOptionType::MustAndOnce, INPUTPARSER::InputType::Parameter);
    // inputParser.AddOption(/* -pa   */ "PROBMAX", "Prob hist max", "100.0", INPUTPARSER::InputNumberRequirement::Equal, 1, INPUTPARSER::InputOptionType::MustAndOnce, INPUTPARSER::InputType::Parameter);
    // inputParser.AddOption(/* -pb   */ "PROBBIN", "Prob hist bin", "0.1", INPUTPARSER::InputNumberRequirement::Equal, 1, INPUTPARSER::InputOptionType::MustAndOnce, INPUTPARSER::InputType::Parameter);
    // inputParser.AddOption(/* -pt   */ "PROBTIME", "Time of prob hist", "0.1", INPUTPARSER::InputNumberRequirement::Equal, 1, INPUTPARSER::InputOptionType::Multiple, INPUTPARSER::InputType::Parameter);
    // inputParser.AddOption(/* -xcr  */ "XCORRELATION", "Coordinate correlation", "false", INPUTPARSER::InputNumberRequirement::Equal, 1, INPUTPARSER::InputOptionType::MustAndOnce, INPUTPARSER::InputType::Parameter);
    // inputParser.AddOption(/* -eta0d*/ "ETA0D", "Eta0 distribution", "false", INPUTPARSER::InputNumberRequirement::Equal, 1, INPUTPARSER::InputOptionType::MustAndOnce, INPUTPARSER::InputType::Parameter);
    // inputParser.AddOption(/* -fh   */ "FPTHIST", "FPT hist min max nbin", "0.0 10.0 100", INPUTPARSER::InputNumberRequirement::Equal, 3, INPUTPARSER::InputOptionType::MustAndOnce, INPUTPARSER::InputType::Parameter);
    if( !inputParser.Parse(&argc, &argv, nodeIndex == 0) ){
        return 0;
    }

    text outputBase = inputParser.GetOptionValue("OUTPUTBASE")[0];
    double x0 = inputParser.GetOptionValue("X0")[0].to_double();
    double v0 = inputParser.GetOptionValue("V0")[0].to_double();
    double gamma = inputParser.GetOptionValue("GAMMA")[0].to_double();
    double kb = inputParser.GetOptionValue("KB")[0].to_double();
    double temperature = inputParser.GetOptionValue("TEMPERATURE")[0].to_double();
    double mass = inputParser.GetOptionValue("MASS")[0].to_double();
    double timestep = inputParser.GetOptionValue("TIMESTEP")[0].to_double();
    size_t numberOfStep = inputParser.GetOptionValue("NUMBEROFSTEPS")[0].to_unsigned_long();
    size_t saveFreq = inputParser.GetOptionValue("FREQSAVE")[0].to_unsigned_long();
    size_t statFreq = inputParser.GetOptionValue("FREQSTAT")[0].to_unsigned_long();
    bool bVr = inputParser.GetOptionValue("RANDOMIZEV")[0].to_boolean();
    bool bEr = true; // inputParser.GetOptionValue("-er")[0].to_boolean();
    bool bSr = inputParser.GetOptionValue("RANDOMIZESTATES")[0].to_boolean();
    bool bSS = inputParser.GetOptionValue("UNIFORMINITALSTATES")[0].to_boolean();
    double periodicLength = inputParser.GetOptionValue("PERIODICLENGTH")[0].to_double();
    double barrierHeight[2][2];
    barrierHeight[0][0] = inputParser.GetOptionValue("BARRIERHEIGHTS0")[0].to_double();
    barrierHeight[0][1] = inputParser.GetOptionValue("BARRIERHEIGHTS0")[1].to_double();
    barrierHeight[1][0] = inputParser.GetOptionValue("BARRIERHEIGHTS1")[0].to_double();
    barrierHeight[1][1] = inputParser.GetOptionValue("BARRIERHEIGHTS1")[1].to_double();
    double initialState = inputParser.GetOptionValue("INITIALSTATE")[0].to_double();
    double slopeHalfLength0 = inputParser.GetOptionValue("BARRIERHALFWIDTH0")[0].to_double();
    double slopeHalfLength1 = inputParser.GetOptionValue("BARRIERHALFWIDTH1")[0].to_double();
    double barrierPosition0 = inputParser.GetOptionValue("BARRIERPOSITION0")[0].to_double();
    double barrierPosition1 = inputParser.GetOptionValue("BARRIERPOSITION1")[0].to_double();
    double q = 0.2; // inputParser.GetOptionValue("-q")[0].to_double();
    double tau = inputParser.GetOptionValue("TAU")[0].to_double();
    double eta0 = 0.0; // inputParser.GetOptionValue("-eta0")[0].to_double();
    bool bWaitingTime = inputParser.GetOptionValue("WAITINGTIME")[0].to_boolean();
    bool bFPT = inputParser.GetOptionValue("FPT")[0].to_boolean();
    bool bTransitionPathTime = inputParser.GetOptionValue("TPT")[0].to_boolean();
    bool bFPTCG = false;
    vector<int> bFPTCGSize;
    if( inputParser.GetOptionSpecified("FPTCG") ){
        for(int j=0;j<inputParser.GetNumberOfOptionValue("FPTCG");j++){
            bFPTCG = inputParser.GetOptionValue("FPTCG", j)[0].to_boolean();
            int value = inputParser.GetOptionValue("FPTCG", j)[1].to_int();
            bFPTCGSize.push_back(value);
        }
    }
    bool bReadRestart = inputParser.GetOptionValue("READRESTART")[0].to_boolean();
    bool bWriteRestart = inputParser.GetOptionValue("WRITERESTART")[0].to_boolean();
    text sInputRestart = inputParser.GetOptionValue("INPUTRESTART")[0];
    text sOutputRestart = inputParser.GetOptionValue("OUTPUTRESTART")[0];
    int ntrial = inputParser.GetOptionValue("NUMBEROFRUNS")[0].to_int();
    size_t numberOfStepPrevious = inputParser.GetOptionValue("NUMBEROFPREVIOUSSTEPS")[0].to_unsigned_long();
    bool bWriteEta = false; // inputParser.GetOptionValue("-we")[0].to_boolean();
    double probHMin = -100.0; // inputParser.GetOptionValue("-pm")[0].to_double();
    double probHMax = 100.0; // inputParser.GetOptionValue("-pa")[0].to_double();
    double probHBin = 0.1; // inputParser.GetOptionValue("-pb")[0].to_double();
    std::vector<double> probTime;
    // for(int i=0;i<inputParser.GetNumberOfOptionValue("-pt");i++){
    //     double value = inputParser.GetOptionValue("-pt", i)[0].to_double();
    //     probTime.push_back(value);
    // }
    bool bWriteX = inputParser.GetOptionValue("WRITEX")[0].to_boolean();
    size_t samplePerOrder = inputParser.GetOptionValue("NUMBEROFLOGSAMPLES")[0].to_unsigned_long();
    int randomNumberSeed = inputParser.GetOptionValue("RANDOMSEED")[0].to_int();
    bool bXCorr = false; // inputParser.GetOptionValue("-xcr")[0].to_boolean();
    bool bTimer = inputParser.GetOptionValue("PRINTTIMER")[0].to_boolean();
    bool bV0D = inputParser.GetOptionValue("DISTRIBUTIONV0")[0].to_boolean();
    bool bEta0D = false; // inputParser.GetOptionValue("-eta0d")[0].to_boolean();
    bool bS0D = inputParser.GetOptionValue("DISTRIBUTIONS0")[0].to_boolean();
    bool bOverdamped = inputParser.GetOptionValue("OVERDAMPED")[0].to_boolean();
    double histBarrierHMin = 0.0; // inputParser.GetOptionValue("-fh")[0].to_double();
    double histBarrierHMax = 10.0; // inputParser.GetOptionValue("-fh")[1].to_double();
    int histBarrierNBin = 100.0; // inputParser.GetOptionValue("-fh")[2].to_int();
    bool bEqualProb = inputParser.GetOptionValue("EQ")[0].to_boolean();
    double equalibriumConstant = inputParser.GetOptionValue("EQ")[1].to_double();


    LangevinEquationRunnerMPI runner(MPI_COMM_WORLD);

    runner.Setup(
        outputBase,
        x0,
        v0,
        gamma,
        kb,
        temperature,
        mass,
        timestep,
        numberOfStep,
        statFreq,
        saveFreq,
        bVr,
        bEr,
        bSr,
        bSS,
        periodicLength,
        barrierHeight[0][0],
        barrierHeight[0][1],
        barrierHeight[1][0],
        barrierHeight[1][1],
        initialState,
        slopeHalfLength0,
        slopeHalfLength1,
        barrierPosition0,
        barrierPosition1,
        q,
        tau,
        eta0,
        bWaitingTime,
        bFPT,
        bReadRestart, 
        bWriteRestart,
        sInputRestart,
        sOutputRestart,        
        ntrial,
        numberOfStepPrevious,
        bWriteEta,
        probHMin,
        probHMax,
        probHBin,
        probTime,
        bWriteX,
        samplePerOrder,
        randomNumberSeed,
        bXCorr,
        bTimer,
        bV0D,
        bEta0D,
        bS0D,
        bOverdamped,
        histBarrierHMin,
        histBarrierHMax,
        histBarrierNBin,
        bTransitionPathTime,
        bFPTCG,
        bFPTCGSize,
        bEqualProb,
        equalibriumConstant);

    runner.Run();

    MPI_Finalize();

}


