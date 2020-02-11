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
  bool null();
  std::string get_name();
  std::vector<std::string> get_route() const;
  Path parent_dir();
  std::string as_str() const;
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
};

namespace IO {
  class Files {
  public:
    Files();
    void initialize(char* argv[]);
    void add_file(std::string name, std::string path, IOtype);

    std::ifstream get_ifstream(std::string name);
    std::ofstream get_ofstream(std::string name);

    bool check_file(std::string name);
    std::string get_file_info(std::string name);

    void print();
    void close();
  private:
    int total_files;
    std::map<std::string, int> file_to_index;
    std::map<int, fileInfo> file_values;
    std::map<int, std::ofstream> ofstream_map;


    Path reference_dir; // Absolute path to location of toml file. Files will be searched for relative to this directory.
    //std::string tomlfile; // Absolute path to options file.
    Path relative_outdir; // Relative path to out dir from reference directory.
    Path absolute_outdir; // Absolute path to the output directory.

    inline std::string path_to_file(int);
    inline void check_stream(const std::string&, const std::string&, std::ifstream&);
    inline std::string findFullFilePath(std::string parameter);
  };
}
#endif
