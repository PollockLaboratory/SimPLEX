#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include <sys/time.h>

using namespace std;




int main() {
	time_t t = time(0);

	int h = t / 3600;
	int m = (t % 3600) / 60;
	int s = t - (t / 60) * 60;

	stringstream time;
	time << setfill('0') << setw(2) << h;
	time << ":" << setfill('0') << setw(2) << m;
	time << ":" << setfill('0') << setw(2) << s;

	cout << time.str() << endl;

	ostringstream ss;
	ss << t;
	// cout << setfill('0') << setw(5) << ss.str();
	cout << setfill('0') << setw(5) << ss.str();


	return (0);
}
