#ifndef __PROGTEST__
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <climits>
#include <cfloat>
#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <string>
#include <vector>
#include <array>
#include <iterator>
#include <set>
#include <list>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <stack>
#include <deque>
#include <memory>
#include <functional>
#include <thread>
#include <mutex>
#include <atomic>
#include <chrono>
#include <stdexcept>
#include <condition_variable>
#include <pthread.h>
#include <semaphore.h>
#include "progtest_solver.h"
#include "sample_tester.h"
using namespace std;
#endif /* __PROGTEST__ */ 

//-------------------------------------------------------------------------------------------------------------------------------------------------------------

class CQualityControl
{
  public:
    CQualityControl();
    static void                        checkAlgorithm                          ( ASheet                                sheet );
    void                               addLine                                 ( AProductionLine                       line );
    void                               start                                   ( int                                   workThreads );
    void                               stop                                    ( void );
  private:
    vector<thread> m_threads;
    int m_num_workThreads;
    int m_CQC_line_id;
    mutex prod_mtx;
    condition_variable prod_cv_empty, cons_cv_empty; 
    struct SSheet
    {
      int sheet_line_id;
      int sheet_id;
      ASheet m_sheet;
    };

    struct SLine
    {
      int line_id;
      mutex line_mtx;
      AProductionLine line;
      condition_variable cv_line_empty;
      atomic_int line_in_sh_id;
      atomic_int line_out_sh_id;
      atomic_bool end;
      atomic_bool fin;
      deque<SSheet> done;
    };

    deque<SLine*> m_lines;

    void in(SLine *& line);
    void out(SLine *& line);

    deque<SSheet> prod_buff;  //comm_in
    void work();

    atomic_bool m_stop;
    atomic_int m_diff;
    // int cons;
};

// TODO: CQualityControl implementation goes here

CQualityControl::CQualityControl()
{
  m_num_workThreads=0;
  m_CQC_line_id=0;
  prod_buff.clear();
  // cons_buff.clear();
  m_lines.clear();
  m_stop=false;
  m_diff=0;
  // cons=0;
}

/*
  the rolling mills are created and registered (method addLine)
*/
/*
  method addLine (x), the method adds a new instance of the rolling mill,
*/
void CQualityControl::addLine ( AProductionLine line )
{
  SLine* a_line = new SLine;
  a_line->line_id=m_CQC_line_id++;
  a_line->line_in_sh_id=0;
  a_line->line_out_sh_id=0;
  a_line->line=line;
  a_line->end=false;
  a_line->fin=false;
  a_line->done.clear();
  m_lines.push_back(a_line);

}

// three conditions to bring get sheet 
// get sheet must not previously return null
// the in_sh_id - out_sh_id < num_workthreads
// *the prod buffer should not be full;
void CQualityControl::in(SLine *& a_line)
{
  while(1)
  {
    // this_thread::sleep_for(1000ms);
    SSheet n_sheet {-1, -1, NULL};
    n_sheet.m_sheet = a_line->line->getSheet();
    ++m_diff;
    if (n_sheet.m_sheet==NULL)
    {
      --m_diff;
      a_line->end=true;
      a_line->cv_line_empty.notify_all();
      break;
    }
    else
    {
      n_sheet.sheet_line_id=a_line->line_id;
      n_sheet.sheet_id=a_line->line_in_sh_id;
      a_line->line_in_sh_id=a_line->line_in_sh_id+1;
      unique_lock<mutex> prod_buff_lock (prod_mtx);
      // DEBUG
      // cout << "in -> prod_buff " << "lineID: "<<n_sheet.sheet_line_id << " sheetID: " << n_sheet.sheet_id << endl;
      // ENDDEBUG
      prod_buff.push_back(n_sheet);
      prod_buff_lock.unlock();
      prod_cv_empty.notify_one();
      // cout << "prod_buff size: " << prod_buff.size() << endl;
    }
  }
}

void CQualityControl::work()
{
  while (1)
  {
    unique_lock<mutex> work_prod_buff_lock (prod_mtx);
    prod_cv_empty.wait(work_prod_buff_lock, [this] () 
    {
      if (!prod_buff.empty())  //return false == wait
      {
        return true;
      }
      else if (m_stop && m_diff==0)
      {
        size_t cnt {0};
        for (size_t i=0; i<m_lines.size(); i++)
        {
          if (m_lines.at(i)->end==false)
          {
            return false;
          }
          else
          {
            cnt++;
          }
          if (cnt == m_lines.size())
          {
            return true;
          }
        }
      }
      else
      {
        return false;
      }
      return false;
    });
    if (prod_buff.empty())
    {
      work_prod_buff_lock.unlock();
      break;
    }
    SSheet a_sheet = prod_buff.front();
    --m_diff;
    prod_cv_empty.notify_all();
    // DEBUG
    // cout << "prod_buff -> worker " << "lineID: " << a_sheet.sheet_line_id << " sheetID: " << a_sheet.sheet_id << endl;
    // ENDDEBUG
    prod_buff.pop_front();
    work_prod_buff_lock.unlock();

    int** arry;
    arry = new int*[a_sheet.m_sheet->m_Length];

    for(int i = 0; i < a_sheet.m_sheet->m_Length; i++) {
      arry[i] = new int[a_sheet.m_sheet->m_Width];
    }

    for(int i = 0; i < a_sheet.m_sheet->m_Length; i++) {
      for(int j = 0; j < a_sheet.m_sheet->m_Width; j++) {
            arry[i][j]=a_sheet.m_sheet->m_Thickness.at((i*a_sheet.m_sheet->m_Width)+j);
        }
    }
    // int loop_inc {0};
    for (map<CRange, int>::iterator it = a_sheet.m_sheet->m_MinMax.begin(); it != a_sheet.m_sheet->m_MinMax.end(); ++it) 
    {
      // int lo= (*it).first.m_Lo;
      // int hi= (*it).first.m_Hi;
      // cout << (*it).first.m_Lo << " " << (*it).first.m_Hi << endl;
      // loop_inc++;
      // int res{-1};
      /*res =*/ (*it).second = maxRectByMinMax(arry, a_sheet.m_sheet->m_Width, a_sheet.m_sheet->m_Length, (*it).first.m_Lo, (*it).first.m_Hi);
      // cout << endl << "sheetID: "<<a_sheet.sheet_id << " " << loop_inc <<" MinMax "<< res << endl;
    }
    // loop_inc=0;
    for (map<int64_t, int>::iterator it = a_sheet.m_sheet->m_Volume.begin(); it != a_sheet.m_sheet->m_Volume.end(); ++it) 
    {
      // int vo= (*it).first;
      // cout << (*it).first << endl;
      // loop_inc++;
      // int res{-1};
      /*res =*/ (*it).second = maxRectByVolume(arry, a_sheet.m_sheet->m_Width, a_sheet.m_sheet->m_Length, (*it).first);
      // cout << endl << "sheetID: "<<a_sheet.sheet_id << " "<< loop_inc <<" Volume "<< res << endl;
    }
    // loop_inc=0;
    for (__cxx11::list<pair<double, int>>::iterator it = a_sheet.m_sheet->m_RelDev.begin(); it != a_sheet.m_sheet->m_RelDev.end(); ++it) 
    {
      // double dev= (*it).first;
      // cout << (*it).first << endl;
      // loop_inc++;
      // int res{-1};
      /*res =*/ (*it).second = maxRectByRelDev(arry, a_sheet.m_sheet->m_Width, a_sheet.m_sheet->m_Length, (*it).first);
      // cout << endl <<"sheetID: "<<a_sheet.sheet_id << " " <<loop_inc <<" RelDev "<< res << endl;
    }

    for(int i = 0; i < a_sheet.m_sheet->m_Length;i++) {
      delete[] arry[i];
    }
    delete[] arry;
    // unique_lock<mutex> work_cons_mutex (cons_mtx);
    // DEBUG
    // cout << "worker -> cons_buff " << "lineID: " << a_sheet.sheet_line_id << " sheetID: " << a_sheet.sheet_id <<endl;
    // ENDDEBUG
    // cons_buff.push_back(a_sheet);
    for (size_t i=0; i<m_lines.size(); i++)
    {
      if (m_lines.at(i)->line_id==a_sheet.sheet_line_id)
      {
        unique_lock<mutex> out_buff_lock (m_lines.at(i)->line_mtx);
        m_lines.at(i)->done.push_back(a_sheet);
        out_buff_lock.unlock();
        m_lines.at(i)->cv_line_empty.notify_one();
        break;
      }
    }
    // work_cons_mutex.unlock();
    // cons++;// TO DELETE
    
    // if ( (m_stop==true) )
    // {
    //   unique_lock<mutex> work_end_lock (prod_mtx);
    //   if (prod_buff.empty())
    //   {
    //     work_end_lock.unlock();
    //     break;
    //   }
    // }
  }
}

void CQualityControl::out(SLine *& a_line)
{
  deque<SSheet> out_c;
  while(1)
  {
    unique_lock<mutex> out_cons_buff_lock (a_line->line_mtx);
    a_line->cv_line_empty.wait(out_cons_buff_lock, [this, &a_line] () 
    {
      if (!(a_line->done.empty()))
      {
        return true;
      }
      else if (a_line->end && a_line->line_in_sh_id==a_line->line_out_sh_id)
      {
        return true;
      }
      else
      {
        return false;
      }
      return false;
    });
    if ( (a_line->end==true) && (a_line->line_out_sh_id==a_line->line_in_sh_id) && (out_c.size()==0))
    {
      out_cons_buff_lock.unlock();
      a_line->fin=true;
      prod_cv_empty.notify_all();
      break;
    }
    out_c.push_back(a_line->done.front());
    // cout << "cons -> out_buff lineID: " << a_line->line_id << "lineID: " << a_line->done.front().sheet_line_id << "sheetID: " << a_line->done.front().sheet_id <<endl;
    a_line->done.pop_front();
    out_cons_buff_lock.unlock();
    sort(out_c.begin(), out_c.end(), [] (SSheet a, SSheet b){return (a.sheet_id<b.sheet_id);});
    while (out_c.size()!=0 && out_c.front().sheet_id==a_line->line_out_sh_id)
    {
      cout << "out_buff -> line;" << a_line->line_id << "line ID: " << out_c.front().sheet_line_id << " sheetID: " << out_c.front().sheet_id << endl ;
      a_line->line->doneSheet(out_c.front().m_sheet);
      a_line->line_out_sh_id++;
      out_c.pop_front();
    }
    
  }
}
/*
  method start ( workThr ), the method starts the communication and worker threads. Once the
  threads are started, method start returns to the caller,
*/
/*
  the computation is started (method start). The method is given the number of worker thread in its parameter.
  Method CQualityControl::start runs the worker threads, and lets them wait for the work. Next, it runs the
  communication threads (two communication threads per registered rolling mill) and lets them serve the mills.
  Once all threads are initialized, the method returns to its caller,
*/
void CQualityControl::start ( int workThreads )
{
  // checkAlgorithm(m_rolling_mill->getSheet());

  m_num_workThreads=workThreads;
  for (int i=0; i<m_num_workThreads; i++)
  {
    m_threads.push_back(thread (&CQualityControl::work, this));
  }
  
  for (long unsigned int i=0; i<m_lines.size(); i++)
  {
    m_threads.push_back(thread(&CQualityControl::in, this, ref(m_lines.at(i))));
    m_threads.push_back(thread (&CQualityControl::out, this, ref(m_lines.at(i))));
  }
}

void CQualityControl::stop ( void )
{
  m_stop=true;
  for (size_t i=0; i<m_threads.size(); i++)
  {
    if ((m_threads.at(i)).joinable())
    {
      (m_threads.at(i)).join();
    }
  }
  while(m_lines.size()!=0)
  {
    delete m_lines.front();
    m_lines.pop_front();
  }
}

void CQualityControl::checkAlgorithm ( ASheet sheet )
{
  // cout << "length: " << sheet->m_Length << ", width: " << sheet->m_Width << ", thickness[0][0]: " << sheet->m_Thickness.at(0) << ", thickness[0][1]: " << sheet->m_Thickness.at(1) << endl;
  // int w = sheet->m_Width;
  // int l = sheet->m_Length;
  // int res {0};
  int** arry;
  arry = new int*[sheet->m_Length];

  for(int i = 0; i < sheet->m_Length; i++) {
    arry[i] = new int[sheet->m_Width];
  }

  for(int i = 0; i < sheet->m_Length; i++) {
    for(int j = 0; j < sheet->m_Width; j++) {
          arry[i][j]=sheet->m_Thickness.at((i*sheet->m_Width)+j);
       }
  }

  for (map<CRange, int>::iterator it = sheet->m_MinMax.begin(); it != sheet->m_MinMax.end(); ++it) 
  {
    // int lo= (*it).first.m_Lo;
    // int hi= (*it).first.m_Hi;
    // cout << (*it).first.m_Lo << " " << (*it).first.m_Hi << endl;
    /*res =*/ (*it).second = maxRectByMinMax(arry, sheet->m_Width, sheet->m_Length, (*it).first.m_Lo, (*it).first.m_Hi);
    // cout << res << endl;
  }

  for (map<int64_t, int>::iterator it = sheet->m_Volume.begin(); it != sheet->m_Volume.end(); ++it) 
  {
    // int vo= (*it).first;
    // cout << (*it).first << endl;
    /*res =*/ (*it).second = maxRectByVolume(arry, sheet->m_Width, sheet->m_Length, (*it).first);
    // cout << res << endl;
  }

  for (__cxx11::list<pair<double, int>>::iterator it = sheet->m_RelDev.begin(); it != sheet->m_RelDev.end(); ++it) 
  {
    // double dev= (*it).first;
    // cout << (*it).first << endl;
    /*res =*/ (*it).second = maxRectByRelDev(arry, sheet->m_Width, sheet->m_Length, (*it).first);
    // cout << res << endl;
  }

  for(int i = 0; i < sheet->m_Length;i++) {
    delete[] arry[i];
  }
  delete[] arry;
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------------
#ifndef __PROGTEST__
int                main                                    ( void )
{
  CQualityControl control;
  AProductionLineTest line0 = std::make_shared<CProductionLineTest> ();
  AProductionLineTest line1 = std::make_shared<CProductionLineTest> ();
  AProductionLineTest line2 = std::make_shared<CProductionLineTest> ();
  control . addLine ( line0 );
  control . addLine ( line1 );
  control . addLine ( line2 );
  control . start ( 2 );
  control . stop  ();
  if ( ! line0 -> allProcessed () )
    throw std::logic_error ( "(some) sheets were not correctly processsed" );
  if ( ! line1 -> allProcessed () )
    throw std::logic_error ( "(some) sheets were not correctly processsed" );
  if ( ! line2 -> allProcessed () )
    throw std::logic_error ( "(some) sheets were not correctly processsed" );
  return 0;  
}
#endif /* __PROGTEST__ */ 
