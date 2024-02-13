// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "NullQMEngine.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "qmmm.h"

#include "SireError/errors.h"

#include "qmmm.h"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

#include "SireError/errors.h"

#include "qmmm.h"

#include "SireError/errors.h"

#include "qmmm.h"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_NullQMEngine_class(){

    { //::SireOpenMM::NullQMEngine
        typedef bp::class_< SireOpenMM::NullQMEngine, bp::bases< SireOpenMM::QMEngine >, boost::noncopyable > NullQMEngine_exposer_t;
        NullQMEngine_exposer_t NullQMEngine_exposer = NullQMEngine_exposer_t( "NullQMEngine", "" );
        bp::scope NullQMEngine_scope( NullQMEngine_exposer );
        { //::SireOpenMM::NullQMEngine::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireOpenMM::NullQMEngine::typeName );
            
            NullQMEngine_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "Get the name of the QM engine." );
        
        }
        { //::SireOpenMM::NullQMEngine::what
        
            typedef char const * ( ::SireOpenMM::NullQMEngine::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireOpenMM::NullQMEngine::what );
            
            NullQMEngine_exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "Get the name of the QM engine." );
        
        }
        NullQMEngine_exposer.staticmethod( "typeName" );
        NullQMEngine_exposer.def( "__str__", &__str__< ::SireOpenMM::NullQMEngine > );
        NullQMEngine_exposer.def( "__repr__", &__str__< ::SireOpenMM::NullQMEngine > );
        NullQMEngine_exposer.def( "__str__", &__str__< ::SireOpenMM::NullQMEngine > );
        NullQMEngine_exposer.def( "__repr__", &__str__< ::SireOpenMM::NullQMEngine > );
    }

}
