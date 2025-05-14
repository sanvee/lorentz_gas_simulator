#pragma once

#include <cmath>

#define TIME_type double

/* ---not needed----
class Time
{
  public:
    double intpart_;
    double fracpart_;

    Time()
    {
      intpart_ = 0.0;
      fracpart_ = 0.0;
    }

    explicit Time ( double x )
    {
      fracpart_ = modf ( x, &intpart_ );
    }

    inline void add ( double x )
    {
      double intpart;
      fracpart_ += modf ( x, &intpart );
      intpart_ += intpart;
    }
    
    inline void add ( Time rhs )
    {
      fracpart_ += rhs.fracpart_;
      intpart_ += rhs.intpart_;
      resplit();
    }
    
    Time (const Time &) = delete;

    inline void resplit ( void )
    {
      if ( intpart_ > 0.0 )
      {
        if ( fracpart_ > 0.0 )
        {
          double temp;
          fracpart_ = modf ( fracpart_, &temp );
          intpart_ += temp;
          return;

        }
        if ( fracpart_ < 0.0 )
        {
          double temp;
          fracpart_ = 1.0 - modf ( fracpart_, &temp );
          intpart_ -= temp + 1.0;
          return;
        }
      }
      else
      {
        if ( fracpart_ < 0.0 )
        {
          double temp;
          fracpart_ = modf ( fracpart_, &temp );
          intpart_ += temp;
          return;
          
        }
        if ( fracpart_ > 0.0 )
        {
          double temp;
          fracpart_ = 1.0 - modf ( fracpart_, &temp );
          intpart_ += temp + 1.0;
          return;
        }
      }
    }

    //copy assignements
    inline Time& operator = ( const Time &rhs )
    {
      fracpart_ = rhs.fracpart_;
      intpart_ = rhs.intpart_;
      return *this;
    }
    inline Time& operator = ( const double rhs )
    {
      fracpart_ = modf ( rhs, &intpart_ );
      return *this;
    }

    inline Time& operator + ( const Time& rhs )
    {
      fracpart_ += rhs.fracpart_;
      intpart_ += rhs.intpart_;
      resplit();
      return *this;
    }

    inline Time& operator - ( const Time& rhs )
    {
      fracpart_ -= rhs.fracpart_;
      intpart_ -= rhs.intpart_;
      resplit();
      return *this;
    }

    inline operator double() const
    {
      return fracpart_ + intpart_;
    }
};

inline bool operator < ( const Time &lhs, const Time &rhs )
{
  if ( lhs.intpart_ == rhs.intpart_ )
  {
    return lhs.fracpart_ < rhs.fracpart_;
  }
  else
  {
    return lhs.intpart_ < rhs.intpart_;
  }
}

inline bool operator> ( const Time& lhs, const Time& rhs )
{
  return rhs < lhs;
}

inline bool operator<= ( const Time& lhs, const Time& rhs )
{
  return ! ( lhs > rhs );
}

inline bool operator>= ( const Time& lhs, const Time& rhs )
{
  return ! ( lhs < rhs );
}

bool operator == ( const Time &lhs, const Time &rhs )
{
  if ( lhs.intpart_ == rhs.intpart_ )
  {
    return lhs.fracpart_ == rhs.fracpart_;
  }
  else
  {
    return lhs.intpart_ == rhs.intpart_;
  }
}

inline bool operator!= ( const Time& lhs, const Time& rhs )
{
  return ! ( lhs == rhs );
}
*/
