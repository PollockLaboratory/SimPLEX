// Input/Output
// Manages reading the input files and writing output files.

#include "Files.h"
#include "../Environment.h"
#include "cpptoml/cpptoml.h"

#include <iostream>
#include <string.h>

#include <sys/stat.h> // For making directories in Linux and OS X
#include <sys/times.h> // For time()
#include <unistd.h> // For getcwd.

extern Environment env;

std::vector<std::string> parse_path(std::string path) {
  // This should error check that the path is correctly formatted.
  std::string rem = path;
  std::vector<std::string> dirs = {};
  while(rem.size() > 0) {
    size_t loc = rem.find("/");

    if(loc == std::string::npos) {
      dirs.push_back(rem);
      break;
    }

    std::string dir = rem.substr(0, loc);
    if(dir != "") {
      dirs.push_back(dir);
    }
    rem = rem.substr(loc + 1, rem.size());
  }
  return(dirs);
}

Path::Path() {
  route = {};
  nullp = false;
}

Path::Path(std::string path) {
  if(path == "") {
    nullp = true;
    route = {};
  } else {
    route = parse_path(path);
    nullp = false;
  }
}

Path::Path(std::vector<std::string> r) {
  if(r.empty() == true) {
    nullp = true;
    route = {};
  } else {
    nullp = false;
    route = r;
  }
}

bool Path::null() {
  return(nullp);
}

std::string Path::get_name() {
  return(route.back());
}

std::vector<std::string> Path::get_route() const {
  return(route);
}

Path Path::parent_dir() {
  // Finds the directory path of the given path to a file.
  std::vector<std::string> new_route = route;
  new_route.pop_back();
  return(Path(new_route));
}

std::string Path::as_str() const {
  auto start = route.begin();
  std::string path = "";
  if(*start == "." or *start == "..") {
    path += *start;
    ++start;
  }
  for(auto it = start; it != route.end(); ++it) {
    path += "/" + *it;
  }
  return(path);
}

std::ostream& operator<<(std::ostream& os, const Path& p) {
  os << p.as_str();
  return(os);
}

Path Path::operator+(const Path& p) {
  if(nullp == true) {
    return(Path(p.get_route()));
  }

  std::vector<std::string> new_route = this->route;
  std::vector<std::string> appending_route = p.get_route();

  for(auto it = appending_route.begin(); it != appending_route.end(); ++it) {
    if(*it == "..") {
      new_route.pop_back();
    } else if(*it != ".") {
      new_route.push_back(*it);
    }
  }

  return(Path(new_route));
}

// Files

IO::Files::Files() {
  total_files = 0;
}

std::string get_current_directory() {
  char cwd[512];
  std::string dir;
  if(getcwd(cwd, sizeof(cwd)) != NULL) {
    dir = std::string(cwd);
  } else {
    std::cerr << "Error: unable to determine current directory." << std::endl;
    exit(EXIT_FAILURE);
  }
  return(dir);
}

void configure_directory(Path dir_path) {
  if(dir_path.null() == true) {
    std::cerr << "Error: unable to configure null directory." << std::endl;
    exit(EXIT_FAILURE);
  }
  
  struct stat st;
  if (stat(dir_path.as_str().c_str(), &st) == 0 && S_ISDIR(st.st_mode)) {
    if (not env.debug) {
      std::cout << "output dir " << dir_path << " exists. Overwrite? (Y/n)" << std::endl;
      if (getchar() != 'Y') exit(1);
    }
  } else {
    mkdir(dir_path.as_str().c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
  }
}

void IO::Files::initialize(char* argv[]) {
  // Set Options file.
  // Find the absolute path to the options file.
  Path cur_dir(get_current_directory());
  Path options_path(argv[1]);

  reference_dir = (cur_dir + options_path).parent_dir();

  if(reference_dir.null() == true) {
    std::cerr << "Error: error specifying reference directory." << std::endl;
    exit(EXIT_FAILURE);
  }
  
  std::cout << "Options file: " << reference_dir + options_path << std::endl;
  add_file("options_file", options_path.get_name(), IOtype::INPUT);

  std::cout << "Command line specified options file used: " << reference_dir + options_path << std::endl;

  relative_outdir = Path(env.get<std::string>("OUTPUT.output_directory"));
  absolute_outdir = reference_dir + relative_outdir;

  std::cout << "Configuring output directory: " << absolute_outdir << std::endl;
  configure_directory(absolute_outdir);
}

inline std::string IO::Files::path_to_file(int i) {
  if(file_values[i].path.null() == true) {
    std::cerr << "Error: trying to get path to NullPath." << std::endl;
    exit(EXIT_FAILURE);
  }
  return((reference_dir + file_values[i].path).as_str());
}

void IO::Files::add_file(std::string name, std::string path, IOtype t) {
  file_to_index[name] = total_files;
  Path fp(path);
  switch(t) {
  case IOtype::INPUT : {

    if(fp.null() == true) {
      std::cerr << "Error: location of required input file not specified." << std::endl;
      exit(EXIT_FAILURE);
    }

    //std::cout << "Input: " << name << " " << fp << std::endl;
    file_values[total_files] = {fp, fp.get_name(), t};
    // Check file exists.
    std::string path = path_to_file(total_files);
    std::ifstream file_stream(path);
    check_stream(name, path, file_stream);
    file_stream.close();
    break;
  }
  case IOtype::OUTPUT : {
    //std::cout << "Output: " << name << " " << path << " " << relative_outdir << std::endl;
    if(path != "") {
      file_values[total_files] = {relative_outdir + fp, fp.get_name(), t};
      ofstream_map[total_files] = get_ofstream(name);
    } else {
      // If no path is given then do not save data.
      file_values[total_files] = {Path(""), "", t};
    }
    break;
  }
  }
  total_files++;
}

std::ifstream IO::Files::get_ifstream(std::string name) {
  int i = file_to_index[name];

  if(file_values[i].t != IOtype::INPUT) {
    std::cout << "Error: Attempting to read " << file_values[i].file_name << " that is not flagged as IOtype::INPUT." << std::endl;
    exit(EXIT_FAILURE);
  }

  std::string path = path_to_file(i);
  std::ifstream file_stream(path);
  check_stream(name, path, file_stream);

  return(file_stream);
}

std::ofstream IO::Files::get_ofstream(std::string name) {
  int i = file_to_index[name];

  if(file_values[i].t != IOtype::OUTPUT) {
    std::cout << "Error: Attempting to read " << file_values[i].file_name << " that is not flagged as IOtype::OUTPUT." << std::endl;
    exit(EXIT_FAILURE);
  }
  
  std::ofstream file_stream(path_to_file(i));
  return(file_stream);
}

// Not quite sure what these are really doing.
bool IO::Files::check_file(std::string name) {
  // Check if an output file is configured not null.
  int i = file_to_index[name];
  Path p = file_values[i].path;
  return(not p.null());
}

void IO::Files::print() {
  std::cout << std::endl << "Files - relative to directory: " << reference_dir << std::endl;

  for(std::map<std::string, int>::iterator it = file_to_index.begin(); it != file_to_index.end(); ++it) {
    if(file_values[it->second].file_name != "") {
      std::cout << it->first << " : " << file_values[it->second].file_name << " ";
      if(file_values[it->second].path.null() == false) {
	std::cout << Path(".") + file_values[it->second].path << " ";
      } else {
	std::cout << ". ";
      }
    } else {
      // When empty path is specified.
      std::cout << it->first << " : None None ";
    }

    if(file_values[it->second].t == IOtype::INPUT) {
      std::cout << "INPUT" << std::endl;
    } else {
      std::cout << "OUTPUT" << std::endl;
    }
  }
  std::cout << std::endl;
}

void copy_file(const Path &sourcefile, const Path &newfile) {
  std::cout << "Copying: ";
  std::cout << sourcefile << " -> " << newfile << std::endl;
  std::ifstream source(sourcefile.as_str());
  std::ofstream destination(newfile.as_str());
  destination << source.rdbuf();
}

void IO::Files::close() {
  for(std::map<std::string, int>::iterator it = file_to_index.begin(); it != file_to_index.end(); ++it) {
    if(file_values[it->second].t == IOtype::INPUT) {
      Path new_path = absolute_outdir + Path(file_values[it->second].file_name);
      copy_file(reference_dir + file_values[it->second].path, new_path);
    }
  }
  print();
}

inline void IO::Files::check_stream(const std::string &name, const std::string &file_path, std::ifstream &s) {
  if (not s.good()) {
    std::cerr << "Error: Cannot read " << name << " file: \"" << file_path << "\"" << std::endl;
    exit(EXIT_FAILURE);
  }
}

std::string IO::Files::get_file_info(std::string name) {
  int i = file_to_index[name];
  return(file_values[i].file_name);
}
