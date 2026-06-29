#pragma once

class LQPLineResult
{
public:
    /*LQPProcessor::Status*/ int status;
    int apex_idx;
    int seed_idx;
    int start_idx;
    int end_idx;
    int code;
    double phys_length;
    double av_field;
    int closed;
    int start, end;
    int pos_min;

    bool useful;

    double *coords;
    int *indices;
    double *distance;

    void init(double *_coords, int *_indices, double *_distance)
    {
        status = 1;
        apex_idx = 0;
        seed_idx = 0;
        start_idx = 0;
        end_idx = 0;
        code = 0;
        phys_length = 0;
        av_field = 0;
        closed = 0;
        start = 0;
        end = 0;
        pos_min = 0;

        useful = false;

        coords = _coords;
        indices = _indices;
        distance = _distance;
    }
};
