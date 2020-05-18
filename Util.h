#pragma once

#include <stdio.h>  /* defines FILENAME_MAX */
#ifdef _WIN32

//the min/max in windows.h conflicts with numeric_limits<T>:max()
//so we need to undefine them
#define NOMINMAX
#include <Windows.h>
#include <direct.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#define GetCurrentDir getcwd
#endif

#include <string>
#include <cctype>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <vector>
#include <set>
#include <stdlib.h>
#include <string.h>

//make sure the configuration is all configuraiton when adding "Additional include directory"!
//otherwise this may say that the file cannot be found.
//#include <boost/filesystem.hpp>



#define uid(x)  ((x) > 0 ? (x) : (-(x)))

namespace CLPLib
{

#ifdef _WIN32


	static int setenv(const char *name, const char *value, int overwrite)
	{
		int errcode = 0;
		if (!overwrite) {
			size_t envsize = 0;
			errcode = getenv_s(&envsize, NULL, 0, name);
			if (errcode || envsize) return errcode;
		}
		return _putenv_s(name, value);
	}

	//used to measure time for both win and linux
	struct timespec { long tv_sec; long tv_nsec; };		//header part

	static int clock_gettime(int, struct timespec *spec)		//C-file part
	{
		__int64 wintime; GetSystemTimeAsFileTime((FILETIME*)&wintime);
		wintime -= 116444736000000000i64;				//1jan1601 to 1jan1970
		spec->tv_sec = (long)(wintime / 10000000i64);           //seconds
		spec->tv_nsec = (long)(wintime % 10000000i64 * 100);      //nano-seconds
		return 0;
	}
#endif

	class Utils
	{
	public:
		static std::string GetCurrentDirectory()
		{
			char cCurrentPath[FILENAME_MAX];
			if (!GetCurrentDir(cCurrentPath, sizeof(cCurrentPath)))
			{
				return "";
			}
			cCurrentPath[sizeof(cCurrentPath) - 1] = '\0'; /* not really required */

			printf("The current working directory is %s", cCurrentPath);
			return std::string(cCurrentPath);
		}

		static std::string& TrimString(std::string& str)
		{
			str.erase(str.begin(), find_if(str.begin(), str.end(),
				[](char& ch)->bool { return !isspace(ch); }));
			str.erase(find_if(str.rbegin(), str.rend(),
				[](char& ch)->bool { return !isspace(ch); }).base(), str.end());
			return str;
		}

		static bool StringEndsWith(const std::string &str, const std::string &suffix)
		{
			return str.size() >= suffix.size() &&
				str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
		}

		static std::vector<std::string> splitString(const std::string &s, char c)
		{
			std::vector<std::string> v;
			std::string::size_type i = 0;
			std::string::size_type j = s.find(c);
			if (j == std::string::npos)
			{
				v.push_back(s);
				return v;
			}
			while (j != std::string::npos)
			{
				v.push_back(s.substr(i, j - i));
				i = ++j;
				j = s.find(c, j);
				if (j == std::string::npos)
				{
					v.push_back(s.substr(i, s.length()));
				}
			}
			return v;
		}


		static std::vector<std::string> splitStringExcludeEmptyString(const std::string &s, char c)
		{
			std::vector<std::string> v;
			v.reserve(16);
			std::string::size_type i = 0;
			std::string::size_type j = s.find(c);
			if (j == std::string::npos)
			{
				if (s.length() > 0)v.push_back(s);
				return v;
			}
			while (j != std::string::npos)
			{
				if (j > i)
				{
					v.push_back(s.substr(i, j - i));
				}
				i = ++j;
				j = s.find(c, j);
				if (j == std::string::npos)
				{
					auto s1 = s.substr(i, s.length());
					if (s1.length() > 0)
					{
						v.push_back(s1);
					}
				}
			}
			return v;
		}

		static int IndexOfAny(const std::string src, const std::vector<char>& anyOf, int startIndex = 0)
		{
			for (unsigned int i = startIndex; i < src.length(); i++)
			{
				for (unsigned int j = 0; j < anyOf.size(); j++)
				{
					if (src[i] == anyOf[j]) return i;
				}
			}
			return -1;
		}

		//similate c# string TryParse
		template<typename T> static bool tryParse(const std::string& s, T& value)
		{
			try
			{
				std::stringstream ss(s);
				if ((ss >> value).fail() || !(ss >> std::ws).eof())
				{
					throw std::bad_cast();
				}
			}
			catch (...)
			{
				return false;
			}
			return true;
		}

		template<typename T> static T parse(const std::string& s)
		{

			T value;
			bool parse_ok = tryParse(s, value);
			if (!parse_ok)
			{
				std::cerr << "An error occured when parsing a valie.\n";
				exit(1);
			}
			return value;
		}

		static void GetPairIDs(long long key, int& idx1, int& idx2)
		{
			idx1 = ((unsigned long long)key) >> 32;
			idx2 = ((unsigned long long)key << 32) >> 32;
		}

		static bool MatchContainsId(long long key, int id)
		{
			int id1, id2;
			GetPairIDs(key, id1, id2);
			return (abs(id1) == abs(id) || abs(id2) == abs(id));
		}

		//given a pair, and id, what is the other id?
		static int GetMatchPairId(long long key, int id)
		{
			int id1, id2;
			GetPairIDs(key, id1, id2);
			int retId = abs(id1) == abs(id) ? abs(id2) : abs(id1);
			return retId;
		}

		//remove parenthesis from the string, e.g. (25)=>25
		static std::string& TrimParenthesis(std::string& str)
		{
			char chars[] = "( )";
			for (unsigned int i = 0; i < strlen(chars); ++i)
			{
				str.erase(std::remove(str.begin(), str.end(), chars[i]), str.end());
			}
			return str;
		}

		//the src vector contains read id, when a read is reverse complemented
		//this read ids change sign.
		static std::vector<int> ReveseOrientation(std::vector<int>& src)
		{
			std::vector<int> vector1;
			vector1.reserve(src.size());
			for (int i : src)
			{
				vector1.push_back(-i);
			}
			return vector1;
		}

		static std::string reverse_complement(const std::string& src)
		{
			auto lambda = [](const char c) {
				switch (c) {
				case 'A':
				case 'a':
					return 'T';
				case 'G':
				case 'g':
					return 'C';
				case 'C':
				case 'c':
					return 'G';
				case 'T':
				case 't':
					return 'A';
				case 'N':
				case 'n':
					return 'N';
				default:
					std::cerr << "warning: reverse complement symbole: " << c << std::endl;
					return c;
				}
			};
			std::string dst;
			dst.resize(src.length());
			std::transform(src.rbegin(), src.rend(), dst.begin(), lambda);
			return dst;
		}

		/*
		static bool fileExists(const std::string path)
		{
			return boost::filesystem::exists(path);
		}
		*/

		static bool RangeOverlap(std::pair<int, int> &a, std::pair<int, int>&b)
		{
			return a.first <= b.second && b.first <= a.second;
		}

		static bool RangeSpan(std::pair<int, int> &a, std::pair<int, int>&b)
		{
			return a.first <= b.first && a.second >= b.second;
		}


		//an array hash function from 
		//http://stackoverflow.com/questions/20511347/a-good-hash-function-for-a-vector
		static std::size_t hash_vector_int(std::vector<int>& vec)
		{
			std::size_t seed = vec.size();
			for (auto& i : vec) {
				seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
			}
			return seed;
		}

		std::size_t hash_set_int(std::set<int>& data) const
		{
			std::size_t seed = data.size();
			for (auto& i : data) {
				seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
			}
			return seed;
		}

		template <class T>
		static std::vector<std::vector<T>> permutate(std::vector<T> &input)
		{
			std::vector<std::vector<T>> output;
			if (input.size() == 1)
			{
				output.push_back(std::vector < T > {input[0]});
			}
			else if (input.size() == 2)
			{
				output.push_back(std::vector < T > {input[0], input[1]});
				output.push_back(std::vector < T > {input[1], input[0]});
			}
			else
			{
				for (int i = 0; i < input.size(); i++)
				{
					T item = input[i];
					std::vector<T> sublist(input);
					sublist.erase(sublist.begin() + i);
					auto subperms = permutate(sublist);
					for (auto& one_perm : subperms)
					{
						one_perm.insert(one_perm.begin(), item);
						output.push_back(one_perm);
					}
				}
			}
			return output;
		}

	};
}