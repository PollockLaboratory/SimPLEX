// Input/Output
// Manages reading the input files and writing output files.

#include "IO.h"
#include "Environment.h"

#include <iostream>
#include <string.h>

#ifdef _WIN32
#include <dir.h> // For mkdir()
#include <ctime> // For time()
#include <Windows.h> // For CreateDirectory. Can't use because of CopyFile
#else
#include <sys/stat.h> // For making directories in Linux and OS X
#include <sys/times.h> // For time()
#endif

extern Environment env;

IO::Files::Files() {
  total_files = 0;
  defaultfile = "resources/defaults.ctrl"; // where to find default settings
  optionsfile = "resources/options.ctrl"; // where to find optional control settings
}

void IO::Files::setupOutputDirectory() {
  outdir = env.get("output_directory");
  ConfigureOutputDirectory();
}

void IO::Files::set_options_file(char* argv[]) {
  std::string s(argv[1]);
  optionsfile = s;
  std::cout << "Command line specified options file used: " << optionsfile << std::endl;
}

void IO::Files::initialize() {
  add_file("default", defaultfile, IOtype::INPUT);
  add_file("options", optionsfile, IOtype::INPUT);
}

void IO::Files::add_file(std::string name, std::string path, IOtype t) {
  file_to_index[name] = total_files;
  std::string file_name = filename_from_path(path);
  switch(t) {
  case IOtype::INPUT : {
    file_values[total_files] = {path, file_name, t};
    // Check file exists.
    std::ifstream file_stream(path);
    check_stream(name, path, file_stream);
    file_stream.close();
    break;
  }
  case IOtype::OUTPUT : {
    path = findFullFilePath(file_name);
    file_values[total_files] = {path, file_name, t};
    break;
  }
  }
  total_files++;
}

std::ifstream IO::Files::get_ifstream(std::string name) {
  int i = file_to_index[name];
  std::string path = file_values[i].path;

  if(file_values[i].t != IOtype::INPUT) {
    std::cout << "Error: Attempting to read " << file_values[i].file_name << " that is not flagged as IOtype::INPUT." << std::endl;
    exit(EXIT_FAILURE);
  }

  std::ifstream file_stream(path);
  check_stream(name, file_values[i].path, file_stream);

  return(file_stream);
}

std::ofstream IO::Files::get_ofstream(std::string name) {
  int i = file_to_index[name];
  std::string path = file_values[i].path;

  if(file_values[i].t != IOtype::OUTPUT) {
    std::cout << "Error: Attempting to read " << file_values[i].file_name << " that is not flagged as IOtype::OUTPUT." << std::endl;
    exit(EXIT_FAILURE);
  }

  std::ofstream file_stream(path);
  return(file_stream);
}

std::string IO::Files::get_file_path(std::string name) {
  int i = file_to_index[name];
  return(file_values[i].path);
}

void IO::Files::print() {
  std::cout << std::endl << "Files:" << std::endl;

  for(std::map<std::string, int>::iterator it = file_to_index.begin(); it != file_to_index.end(); ++it) {
    std::cout << it->first << " : " << file_values[it->second].file_name << " " << file_values[it->second].path << " ";
    if(file_values[it->second].t == IOtype::INPUT) {
      std::cout << "INPUT" << std::endl;
    } else {
      std::cout << "OUTPUT" << std::endl;
    }
  }
  std::cout << std::endl;
}

void IO::Files::close() {
  for(std::map<std::string, int>::iterator it = file_to_index.begin(); it != file_to_index.end(); ++it) {
    if(file_values[it->second].t == IOtype::INPUT) {
      std::string new_path = findFullFilePath(file_values[it->second].file_name);
      copyFile(file_values[it->second].path, new_path);
    }
  }
  print();
}

void IO::Files::check_stream(const std::string &name, const std::string &file_path, std::ifstream &s) {
  if (not s.good()) {
    std::cerr << "Cannot read " << name << " file: \"" << file_path << "\"" << std::endl;
    exit(EXIT_FAILURE);
  }
}

std::string IO::Files::filename_from_path(std::string path) {
  char * cstr = new char[path.length() + 1];
  strcpy(cstr, path.c_str());

  char *pch;
  pch = strtok(cstr, "/");
  char *lastToken;

  while(pch != NULL) {
    lastToken = pch;
    pch = strtok(NULL, "/");
  }
  return(std::string(lastToken));
}

void IO::Files::copyFile(const string &sourcefile, const string &newfile) {
  std::cout << "Copying" << std::endl;
  std::cout << sourcefile << " " << newfile << std::endl;
  std::ifstream source(sourcefile);
  std::ofstream destination(newfile);
  destination << source.rdbuf();
}

inline std::string IO::Files::findFullFilePath(std::string parameter) {
  /*
   * Prepends the output directory path to the output file names, giving the full path name
   * for the given file.
   */
  if (parameter == "") {
    std::cerr << "Cannot prepend output directory to empty parameter" << std::endl;
  }
  parameter = outdir + parameter;
  return parameter;
}

void IO::Files::ConfigureOutputDirectory() {
  /*
   * Configures the output directory and the file names of the output files.
   * Creates the output directory and changes the file names of the output file to
   * include to full path to the output directory.
   *
   * For example: "seq.out" -> "/output_dir/seq.out"
   */

  char lastchar = outdir.at(outdir.length() - 1);

  if(lastchar != '/' && lastchar != '\\') {
    if(env.debug) std::cout << "last char is not /" << std::endl;
    outdir += '/';
  }

  // For Windows
  // According to http://sourceforge.net/p/predef/wiki/OperatingSystems/
  // TODO: double check this works on windows for directory that already exists
#ifdef _WIN32
  if (mkdir(outdir.c_str())) {
    std::cout << "Could not make output directory " << output_directory << std::endl;
  } else {
    std::cout << "Making directory successful" << std::endl;
  }
#else	//ifdef __linux__ || __APPLE__
  struct stat st;
  if (stat(outdir.c_str(), &st) == 0 && S_ISDIR(st.st_mode)) {
    if (not env.debug) {
      std::cout << "output dir " << outdir << " exists. Overwrite? (Y/n)" << std::endl;
      if (getchar() != 'Y') exit(1);
    }
  } else {
    mkdir(outdir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
  }

  /* Read, write, and search, or execute, for the file owner; S_IRWXU is the bitwise inclusive
   * OR of S_IRUSR, S_IWUSR, and S_IXUSR.
   * Read, write, and search or execute permission for the file's group. S_IRWXG is the
   * bitwise inclusive OR of S_IRGRP, S_IWGRP, and S_IXGRP. Read, write, and search or
   * execute permission for users other than the file owner. S_IRWXO is the bitwise inclusive
   * OR of S_IROTH, S_IWOTH, and S_IXOTH.*/
}
#endif
