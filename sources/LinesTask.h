#pragma once

#include "TaskQueueProcessor.h"

//---------------------------------------------------------------------------------------
class LQPTask : public ATQPTask
{
protected:
    double data[3];
    unsigned long queueID;

public:
    LQPTask(int _id) : ATQPTask(_id) {}
    void setData(double *_data) { memcpy(data, _data, 3 * sizeof(double)); }
    double *getData() { return data; }
    void setID(unsigned long _id) { queueID = _id; }
    unsigned long getID() { return queueID; }
};

//---------------------------------------------------------------------------------------
class LQPTaskFactory : public ATQPTaskFactory
{
public:
    LQPTaskFactory() : ATQPTaskFactory() {}
    virtual ~LQPTaskFactory() {}

    virtual ATQPTask *create()
    {
        return new LQPTask(counter++);
    }
};
