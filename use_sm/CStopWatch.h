// CStopWatch.h
// Created by Fares Alhassen, 2010.08.04
// Description: Provides a simple class to measure elasped computational times in C++ application

#ifndef CSTOPWATCH_H_INCLUDED
#define CSTOPWATCH_H_INCLUDED

#include <time.h>
#include <iostream>

// define structures
// time variable
typedef clock_t TimeVariable;

class CStopWatch {
    public:
        // constructors/destructors
        CStopWatch();
        ~CStopWatch();

        // methods
        void start();
        void measure();
        void reset();

    private:
        // double startTimeInSeconds, elapsedTimeInSeconds;
        TimeVariable startTime, elapsedTime, lastMeasuredTime;
        double getTimeDifferenceInSeconds( TimeVariable time2, TimeVariable time1 );
        double getTimeDifferenceInMinutes( TimeVariable time2, TimeVariable time1 );
};

// constructor
CStopWatch::CStopWatch() {
}

// destructor
CStopWatch::~CStopWatch() {
}

// methods

// starts the stopwatch
void CStopWatch::start() {
    startTime = clock();
    lastMeasuredTime = startTime;
}

// measures the elapsed time on the stopwatch
void CStopWatch::measure() {

    elapsedTime = clock();

    std::cout << "Elapsed time: " << getTimeDifferenceInMinutes( elapsedTime, lastMeasuredTime ) << " min / ";
    std::cout << getTimeDifferenceInMinutes( elapsedTime, startTime ) << " min ( ";
    std::cout << getTimeDifferenceInSeconds( elapsedTime, lastMeasuredTime ) << " sec / ";
    std::cout << getTimeDifferenceInSeconds( elapsedTime, startTime ) << " sec ) " << std::endl << std::endl;

    // storing current elapsed time as last measurement
    lastMeasuredTime = elapsedTime;

}

// resets the stopwatch
void CStopWatch::reset() {
    start();
}

// outputs the time difference
double CStopWatch::getTimeDifferenceInSeconds( TimeVariable time2, TimeVariable time1 ) {
    return( static_cast< double >( ( time2 - time1 ) / static_cast< double >( CLOCKS_PER_SEC ) ) );
}

double CStopWatch::getTimeDifferenceInMinutes( TimeVariable time2, TimeVariable time1 ) {
    return( getTimeDifferenceInSeconds( time2, time1 ) / 60.0 );
}

#endif
