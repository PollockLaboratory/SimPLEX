#ifndef Files_h_
#define Files_h_

#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <iostream>

// Path classes.

class Path {
  // Base class to file or directory.
private:
  bool nullp;
  std::vector<std::string> route;
public:
  Path();
  Path(std::string);
  Path(std::vector<std::string>);
  bool null() const;
  std::string get_name();
  std::vector<std::string> get_route() const;
  Path parent_dir();
  std::string as_str() const;
  void add_timestamp();
  friend std::ostream& operator<<(std::ostream& os, const Path& path);
  Path operator+(const Path&);
};

enum class IOtype {
  INPUT,
  OUTPUT
};

struct fileInfo {
  Path path;
  std::string file_name;
  IOtype t;
  int fd;
  int size;
};

namespace IO {
  class Files {
  public:
    Files();
    void initialize(char* argv[]);

    void add_file(std::string name, std::string path, IOtype);
    fileInfo& get_info(std::string name);

    bool end(std::string name);
    std::string get_next_line(std::string name);
    std::string read_all(std::string name);

    void write_to_file(std::string name, std::string data);
    
    void print();
    void clean_and_close();
  private:
    int total_files;
    std::map<std::string, int> file_to_index;
    std::map<int, fileInfo> file_values;

    Path reference_dir; // Absolute path to location of options toml file.
                        // Files will be searched for relative to this directory.
    Path relative_outdir; // Relative path to out dir from reference directory.
    Path absolute_outdir; // Absolute path to the output directory.
  };
}
#endif
