#include "CRADLE/ConfigParser.hh"

#include "CLI11.hpp"

namespace CRADLE
{

  void parse(CLI::App &app, int argc, const char **argv)
  {
    if (argc == 0)
    {
      argc = 1;
      const char *_argv[1] = {"coupling"};
      argv = _argv;
    }
    try
    {
      app.parse(argc, argv);
    }
    catch (const CLI::ParseError &e)
    {
      app.exit(e);
    }
  }

  void SetGeneralOptions(CLI::App &app, General &general)
  {
    CLI::App *cmd = app.add_subcommand("General", "This is the general subcommand")->ignore_case()->required();
    cmd->add_option("-v,--Verbosity", general.Verbosity, "Verbosity settings");
    cmd->add_option("-f,--Verbosity_file", general.Verbosity_file, "Verbosity File settings");
    cmd->add_option("-l,--loop", general.Loop, "Number of events to generate.")->required();
    cmd->add_option("-t,--threads", general.Threads, "Number of threads (2 x #CPU).");
    cmd->add_option("-o,--output", general.Output, "Name of the output file.");
  }

  void SetNuclearOptions(CLI::App &app, NuclearOptions &nuclear)
  {

    CLI::App *cmd = app.add_subcommand("Nucleus", "This is the initial nucleus options command.")->ignore_case()->required();
    cmd->add_option("-n,--name", nuclear.Name, "Name of initial particle.")->required();
    cmd->add_option("-Z,--charge", nuclear.Charge, "Charge as multiple of proton charge.");
    cmd->add_option("-A,--mass", nuclear.Nucleons, "Mass number.");
    cmd->add_option("-e,--energy", nuclear.Energy, "Excitation energy of initial state.");
    cmd->add_option("-w,--weakmagnetism", nuclear.WeakMagnetism, "Weak magnetism term bA/c.");
    cmd->add_option("--Alignment", nuclear.Alignment, "");
    cmd->add_option("--PolarisationMag", nuclear.PolarisationMag, "");
    cmd->add_option("--PolarisationX", nuclear.PolarisationX, "");
    cmd->add_option("--PolarisationY", nuclear.PolarisationY, "");
    cmd->add_option("--PolarisationZ", nuclear.PolarisationZ, "");
  }

  void SetCouplingConstants(CLI::App &app, CouplingConstants &couplingConstants)
  {
    CLI::App *cmd = app.add_subcommand("Coupling", "This is the coupling subcommand");
    cmd->add_option("--CV", couplingConstants.CV, "Vector coupling constant.");
    cmd->add_option("--CT", couplingConstants.CT, "Tensor coupling constant.");
    cmd->add_option("--CS", couplingConstants.CS, "Scalar coupling constant.");
    cmd->add_option("--CA", couplingConstants.CA, "Axial coupling constant.");
    cmd->add_option("--CVP", couplingConstants.CVP, "Vector prime coupling constant.");
    cmd->add_option("--CTP", couplingConstants.CTP, "Tensor prime coupling constant.");
    cmd->add_option("--CSP", couplingConstants.CSP, "Scalar prime coupling constant.");
    cmd->add_option("--CAP", couplingConstants.CAP, "Axial prime coupling constant.");
    cmd->add_option("--a", couplingConstants.a, "Angular correlation coefficient.");
    cmd->add_option("--b", couplingConstants.b, "Fierz term.");
    cmd->add_option("--A", couplingConstants.A, "Beta asymmetry.");
    cmd->add_option("--B", couplingConstants.B, "Neutrino asymmetry.");
    cmd->add_option("--D", couplingConstants.D, "Triple correlation coefficient.");
    cmd->add_option("--c", couplingConstants.c, "Energy dependence of the correlation coefficients.");
  }

  void SetCuts(CLI::App &app, Cuts &cuts)
  {
    CLI::App *cmd = app.add_subcommand("Cuts", "This is the cuts subcommand")->ignore_case();
    cmd->add_option("--distance", cuts.Distance, "")->ignore_case();
    cmd->add_option("--lifetime", cuts.Lifetime, "")->ignore_case();
    cmd->add_option("--energy", cuts.Energy, "")->ignore_case();
  }

  void SetBetaDecayOptions(CLI::App &app, BetaDecay &betaDecay)
  {
    CLI::App *cmd = app.add_subcommand("BetaDecay", "This is the beta decay subcommand")->ignore_case();
    cmd->add_option("--Default", betaDecay.Default, "");
    cmd->add_option("--FermiFunction", betaDecay.FermiFunction, "");
  }

  void SetDecayOptions(CLI::App &app, Decay &decay)
  {
    CLI::App *cmd = app.add_subcommand("Decay", "This is the decay subcommand")->ignore_case();
    cmd->add_option("--Nuclear_Level_Width", decay.Nuclear_Level_Width, "")->ignore_case();
    cmd->add_option("--InFlightDecay", decay.InFlightDecay, "")->ignore_case();
    cmd->add_option("--GammaGammaCorrelation", decay.GammaGammaCorrelation, "")->ignore_case();
  }

  void SetEnvironmentOptions(CLI::App &app, EnvOptions &envOptions)
  {
    app.add_option("--AMEdata", envOptions.AMEdata, "AME2020 file location")->envname("AMEdata");
    app.add_option("--Gammadata", envOptions.Gammadata, "")->envname("Gammadata");
    app.add_option("--Radiationdata", envOptions.Radiationdata, "")->envname("Radiationdata");
  }

  ConfigOptions ParseOptions(std::string filename, int argc, const char **argv)
  {
    ConfigOptions configOptions;

    CLI::App app{"CRADLE++ App"};
    app.allow_extras(true);
    app.allow_config_extras(true);
    app.set_config("-c", filename);

    SetGeneralOptions(app, configOptions.general);
    SetNuclearOptions(app, configOptions.nuclear);
    SetCouplingConstants(app, configOptions.couplingConstants);
    SetCuts(app, configOptions.cuts);
    SetBetaDecayOptions(app, configOptions.betaDecay);
    SetDecayOptions(app, configOptions.decay);
    SetEnvironmentOptions(app, configOptions.envOptions);

    parse(app, argc, argv);

    PrintingAllOptions(configOptions);

    std::cout << std::endl;

    return configOptions;
  }

  void PrintingAllOptions(const ConfigOptions &configOptions)
  {
    Message("General", "", 0, "CYAN");
    Message("General", "Verbosity: " + std::to_string(configOptions.general.Verbosity), 1, "blue");
    Message("General", "Verbosity file: " + std::to_string(configOptions.general.Verbosity_file), 1, "blue");
    Message("General", "Loop: " + std::to_string(configOptions.general.Loop), 1, "blue");
    Message("General", "Threads: " + std::to_string(configOptions.general.Threads), 1, "blue");
    Message("General", "Output: " + configOptions.general.Output, 1, "blue");

    Message("Nucleus", "", 0, "CYAN");
    Message("Nucleus", "Name: " + configOptions.nuclear.Name, 1, "blue");
    Message("Nucleus", "Z: " + std::to_string(configOptions.nuclear.Charge), 1, "blue");
    Message("Nucleus", "A: " + std::to_string(configOptions.nuclear.Nucleons), 1, "blue");
    Message("Nucleus", Form("Energy: %.1f keV", configOptions.nuclear.Energy), 1, "blue");
    Message("Nucleus", Form("bWM: %.3f", configOptions.nuclear.WeakMagnetism), 1, "blue");
    Message("Nucleus", Form("Alignment: %.2f", configOptions.nuclear.Alignment), 1, "blue");
    Message("Nucleus", Form("PolarisationMag: %.2f", configOptions.nuclear.PolarisationMag), 1, "blue");
    Message("Nucleus", Form("PolarisationX: %.2f", configOptions.nuclear.PolarisationX), 1, "blue");
    Message("Nucleus", Form("PolarisationY: %.2f", configOptions.nuclear.PolarisationY), 1, "blue");
    Message("Nucleus", Form("PolarisationZ: %.2f", configOptions.nuclear.PolarisationZ), 1, "blue");

    Message("Coupling", "", 0, "CYAN");
    Message("Coupling", Form("Cᵥ: (%.2f%+.2fi)", configOptions.couplingConstants.CV.real(), configOptions.couplingConstants.CV.imag()), 1, "blue");
    Message("Coupling", Form("Cᵥ': (%.2f%+.2fi)", configOptions.couplingConstants.CVP.real(), configOptions.couplingConstants.CVP.imag()), 1, "blue");
    Message("Coupling", Form("Cₐ: (%.2f%+.2fi)", configOptions.couplingConstants.CA.real(), configOptions.couplingConstants.CA.imag()), 1, "blue");
    Message("Coupling", Form("Cₐ': (%.2f%+.2fi)", configOptions.couplingConstants.CAP.real(), configOptions.couplingConstants.CAP.imag()), 1, "blue");
    Message("Coupling", Form("Cₛ: (%.2f%+.2fi)", configOptions.couplingConstants.CS.real(), configOptions.couplingConstants.CS.imag()), 1, "blue");
    Message("Coupling", Form("Cₛ': (%.2f%+.2fi)", configOptions.couplingConstants.CSP.real(), configOptions.couplingConstants.CSP.imag()), 1, "blue");
    Message("Coupling", Form("Cₜ: (%.2f%+.2fi)", configOptions.couplingConstants.CT.real(), configOptions.couplingConstants.CT.imag()), 1, "blue");
    Message("Coupling", Form("Cₜ': (%.2f%+.2fi)", configOptions.couplingConstants.CTP.real(), configOptions.couplingConstants.CTP.imag()), 1, "blue");
    Message("Correlation", "", 0, "CYAN");
    Message("Correlation", Form("aᵦᵥ: %.2f", configOptions.couplingConstants.a), 1, "blue");
    Message("Correlation", Form("bᶠ: %.2f", configOptions.couplingConstants.b), 1, "blue");
    Message("Correlation", Form("Aᵦ: %.2f", configOptions.couplingConstants.A), 1, "blue");
    Message("Correlation", Form("Bᵥ: %.2f", configOptions.couplingConstants.B), 1, "blue");
    Message("Correlation", Form("D: %.2f", configOptions.couplingConstants.D), 1, "blue");
    Message("Correlation", Form("c: %.2f", configOptions.couplingConstants.c), 1, "blue");

    Message("Cuts", "", 0, "CYAN");
    Message("Cuts", Form("Distance: %.2f", configOptions.cuts.Distance), 1, "blue");
    Message("Cuts", Form("Lifetime: %.2f", configOptions.cuts.Lifetime), 1, "blue");
    Message("Cuts", Form("Energy: %.2f", configOptions.cuts.Energy), 1, "blue");

    Message("BetaDecay", "", 0, "CYAN");
    Message("BetaDecay", "Default: " + configOptions.betaDecay.Default, 1, "blue");
    Message("BetaDecay", "FermiFunction: " + configOptions.betaDecay.FermiFunction, 1, "blue");

    Message("Decay", "", 0, "CYAN");
    Message("Decay", Form("InFlightDecay: %s", configOptions.decay.InFlightDecay ? "true" : "false"), 1, "blue");
    Message("Decay", Form("Nuclear_Level_Width: %s", configOptions.decay.Nuclear_Level_Width ? "true" : "false"), 1, "blue");
    Message("Decay", Form("γγ Correlation: %s", configOptions.decay.GammaGammaCorrelation ? "true" : "false"), 1, "blue");
  }
} // end of namespace CRADLE
