#pragma once
#include "Exception.h"

#include <chrono>
#include <iostream>
#include <windows.h>	//compile error is occured when windows.h is not includes before psapi.h
#include <psapi.h>
#include <vector>

class Profiler
{
public:
	static void set_time_point(void);
	static double get_time_duration(void);

	static void record_Consumed_Memory(void);
	static void print_Consumed_Memory(void);
	static void record_Consumed_Memory_And_Time(void);
	static void print_Consumed_Memory_And_Time(void);	

private:
	static inline PROCESS_MEMORY_COUNTERS_EX memory_recorder_;
	static inline std::vector<size_t> memory_record_;
	static inline std::vector<std::chrono::steady_clock::time_point> time_points_;
};