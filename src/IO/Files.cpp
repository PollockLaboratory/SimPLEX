// Input/Output
// Manages reading the input files and writing output files.

#include "Files.h"
#include "../Environment.h"
#include "cpptoml/cpptoml.h"

#include <iostream>
#include <ctime>
#include <string.h>

#include <sys/types.h>
#include <sys/stat.h> // For making directories in Linux and OS X
#include <unistd.h>
#include <fcntl.h>
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

bool Path::null() const {
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

void Path::add_timestamp() {
  auto it = --route.end();
  std::time_t t = std::time(0);
  std::tm* now = std::localtime(&t);

  char buffer[50];
  sprintf(buffer, "_%i-%i-%i_%i:%i", (now->tm_year + 1900), (now->tm_mon + 1), now->tm_mday, now->tm_hour, now->tm_min);

  std::cout << "Last element: " << *it << std::endl;
  *it += buffer;
  std::cout << "Last element: " << *it << std::endl;
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

fileInfo& IO::Files::get_info(std::string name) {
  if(file_to_index.find(name) == file_to_index.end()) {
    std::cerr << "Error: \'" << name << "\' is being requested, but has not been added to Files." << std::endl;
    exit(EXIT_FAILURE);
  }
  int id = file_to_index[name];
  return(file_values[id]);
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
  if(stat(dir_path.as_str().c_str(), &st) == 0 && S_ISDIR(st.st_mode)) {
    if(not env.get<bool>("OUTPUT.overwrite_output")) {
      std::cout << "Error: cannot write to the specified output directory, as it already exists." << std::endl;
      std::cout << "Either change the specified output directory or set OUTPUT.overwrite_output to True." << std::endl;
      exit(EXIT_FAILURE);
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

  //std::string outdir_name = env.get<std::string>("OUTPUT.output_directory");

  relative_outdir = Path(env.get<std::string>("OUTPUT.output_directory"));
  if(env.get<bool>("OUTPUT.output_directory_append_time")) {
    relative_outdir.add_timestamp();
  }

  //if(env.get<bool>("OUTPUT.output
  absolute_outdir = reference_dir + relative_outdir;

  std::cout << "Configuring output directory: " << absolute_outdir << std::endl;
  configure_directory(absolute_outdir);
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

    std::string full_path = (reference_dir + fp).as_str();

    int fd = open(full_path.c_str(), O_RDONLY);
    if(fd == -1) {
      std::cerr << "Error: unable to open \'" << path << "\'." << std::endl;
      exit(EXIT_FAILURE);
    }
    struct stat finfo;
    fstat(fd, &finfo);
    int size = finfo.st_size;

    file_values[total_files] = {fp, fp.get_name(), t, fd, size};
    break;
  }
  case IOtype::OUTPUT : {
    //std::cout << "Output: " << name << " " << path << " " << relative_outdir << std::endl;
    if(path != "") {
      std::string full_path = (absolute_outdir + fp).as_str();
      mode_t file_permissions = S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH;
      int fd = open(full_path.c_str(), O_WRONLY | O_CREAT | O_TRUNC, file_permissions);
      file_values[total_files] = {absolute_outdir + fp, fp.get_name(), t, fd, 0};
    } else {
      // If no path is given then do not save data.
      file_values[total_files] = {Path(""), "", t, 0, 0};
    }
    break;
  }
  }
  total_files++;
}

bool IO::Files::end(std::string name) {
  fileInfo info = get_info(name);

  if(lseek(info.fd, 0, SEEK_CUR) >= info.size) {
    return(true);
  } else {
    return(false);
  }
}

std::string IO::Files::get_next_line(std::string name) {
  fileInfo info = get_info(name);

  if(info.t != IOtype::INPUT) {
    std::cerr << "Error: Attempting to read " << info.file_name << " that is not flagged as IOtype::INPUT." << std::endl;
    exit(EXIT_FAILURE);
  }

  std::string data = "";
  ssize_t num_read;

  int fd = info.fd;
  int buf_size = 256;
  char buf[buf_size];

  while(true) {
    int loc = lseek(fd, 0, SEEK_CUR);
    num_read = read(fd, &buf, buf_size);

    if(num_read == -1) {
      std::cerr << "Error: cannot read \'" << info.file_name << "\'." << std::endl;
      exit(EXIT_FAILURE);
    }

    for(int i = 0; i < buf_size; i++) {
      if(loc + i >= info.size) {
	return(data);
      }

      if(buf[i] == '\n') {
	lseek(fd, i-num_read+1, SEEK_CUR);
	return(data);
      } else {
	data.push_back(buf[i]);
      }
    }
  }  
}

std::string IO::Files::read_all(std::string name) {
  fileInfo info = get_info(name);
  lseek(info.fd, 0, SEEK_SET);

  char* buf = new char [info.size];
  ssize_t num_read = read(info.fd, buf, info.size);

  if(num_read == -1) {
    std::cerr << "Error: cannot read \'" << info.file_name << "\'." << std::endl;
    exit(EXIT_FAILURE);
  }

  std::string ret(buf, info.size);

  delete[] buf;

  return(ret);
}

void IO::Files::write_to_file(std::string name, std::string data) {
  fileInfo info = get_info(name);

  if(info.file_name == "") {
    return;
  }

  const char* cstr = data.c_str();
  ssize_t num_write = write(info.fd, cstr, data.size());

  if(num_write == -1) {
    std::cerr << "Error: failed to write to \'" << info.file_name << "\'." << std::endl;
    exit(EXIT_FAILURE);
  }

  // Increase size of file in fileInfo.
}

void IO::Files::print() {
  unsigned int max_name_len = 0;
  unsigned int max_file_len = 4; // If empty no file_name then 'None' is printed, which is len = 4.

  for(std::map<std::string, int>::iterator it = file_to_index.begin(); it != file_to_index.end(); ++it) {
    it->first.length() > max_name_len ? max_name_len = it->first.length() : max_name_len = max_name_len;
    unsigned int file_len = file_values[it->second].file_name.length();
    file_len > max_file_len ? max_file_len = file_len : max_name_len = max_name_len;
  }

  max_name_len++;
  max_file_len++;

  std::cout << std::endl << "Files - relative to directory: " << reference_dir << std::endl;

  for(std::map<std::string, int>::iterator it = file_to_index.begin(); it != file_to_index.end(); ++it) {
    // Reference name.
    std::cout << it->first;
    for(unsigned int i = 0; i < max_name_len - it->first.length(); i++) {
      std::cout << " ";
    }

    if(file_values[it->second].t == IOtype::INPUT) {
      std::cout << "[INPUT ] ";
    } else {
      std::cout << "[OUTPUT] ";
    }

    // File name.
    std::string file_name;
    if(file_values[it->second].file_name != "") {
      file_name = file_values[it->second].file_name;
    } else {
      file_name = "None";
    }

    std::cout << file_name;
    for(unsigned int i = 0; i < max_file_len - file_name.length(); i++) {
      std::cout << " ";
    }

    // File path
    if(file_values[it->second].path.null() == false) {
      std::cout << Path(".") + file_values[it->second].path;
    } else {
	std::cout << "None";
    }

    std::cout << std::endl;
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

void IO::Files::clean_and_close() {
  for(std::map<std::string, int>::iterator it = file_to_index.begin(); it != file_to_index.end(); ++it) {
    if(file_values[it->second].t == IOtype::INPUT) {
      Path new_path = absolute_outdir + Path(file_values[it->second].file_name);
      copy_file(reference_dir + file_values[it->second].path, new_path);
    }
    fileInfo info = get_info(it->first);
    close(info.fd);
  }
  print();
}
