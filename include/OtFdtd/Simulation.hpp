/* Fdtd/Simulation.hpp - The object for running a FDTD simulation.
 *
 * Copyright (C) 2016 Isaac Lenton (aka ilent2)
 */

#ifndef FDTD_SIMULATION_HPP
#define FDTD_SIMULATION_HPP

#include "Modules/EmptyModule.hpp"
#include "Modules/FieldData.hpp"
#include "Modules/ModuleSelectors.hpp"
#include "Utilities/DecoratedChild.hpp"
#include "Utilities/DecoratedLayers.hpp"

namespace Fdtd {

/** An object for running FDTD simulations.
 *
 * All simulations must have Timing and Indices modules, the ElementData
 * members of these modules are ignored.
 *
 * @tparam Timing        : Module describing simulation time evolution.
 * @tparam Indices       : Module describing simulation time evolution.
 * @tparam Decorators... : Layers for adding Permittivyt, TFSF, CPML, etc.
 */
template <typename Timing, typename Indices, typename... Decorators>
class Simulation
{
  //
  // Resolve FieldData type
  //

  /** Base type used for element data. */
  class ElementDataBase {};

  /** Assembled element data type. */
  using ElementData = typename Utilities::DecoratedChild<
      ElementDataSelector, Decorators...>::template Type<ElementDataBase>;

  /** The field data type to be added after Timing and Indices. */
  using FieldData = Fdtd::FieldDataModule<ElementData>;

  //
  // Base type for Field type
  //

  /** Base type used for Field type. */
  struct FieldBase : Fdtd::EmptyModule
  {
    /** Describes how the simulation performed. */
    template <typename Base>
    struct Evaluate : public Base
    {
      /** Default method for running the simulation. */
      void evaluate(void)
      {
        Base::initialize();

        for (auto it = Base::TimingBegin(); it != Base::TimingEnd(); ++it)
        {
          Base::evolve(it);
        }

        Base::finalize();
      }
    };

    /** Fallback methods for initialize/finalize. */
    template <typename Base>
    struct Fallback : public Base
    {
      /** Terminal case for module initialization.
       *
       * Included in Fallback so any module/sub-module can implement.
       *
       * A module implementing this method should first call
       * Base::initialize and then perform its own initialization.
       */
      void initialize(void)
      {
      }

      /** Terminal case for module finalization.
       *
       * Included in Fallback so any module/sub-module can implement.
       *
       * A module implementing this method should perform its own
       * finalization before calling Base::finalize.
       */
      void finalize(void)
      {
      }
    };

    /** Provides the methods used by evaluate. */
    template <typename Base>
    struct EvaluateMethods : public Base
    {
      /** Default 2-step evolution implementation. */
      void evolve(typename Base::TimingIterator it)
      {
        // Evolve H-field and write state
        Base::evolveH(it);
        Base::updateH(it);

        // Evolve E-field and write state
        Base::evolveE(it);
        Base::updateE(it);
      }
    };

    /** Provides the 2-step methods for advancing the simulation. */
    template <typename Base>
    struct Evolve : public Base
    {
      /** Default for evolving H-field with a 2-step method. */
      void evolveH(typename Base::TimingIterator it)
      {
        Base::advanceH(it);
      }

      /** Default for writing post H-field evolution. */
      void updateH(typename Base::TimingIterator it)
      {
        Base::inspectH(it);
      }

      /** Default for evolving E-field with a 2-step method. */
      void evolveE(typename Base::TimingIterator it)
      {
        Base::advanceE(it);
      }

      /** Default for writing post E-field evolution. */
      void updateE(typename Base::TimingIterator it)
      {
        Base::inspectE(it);
      }
    };

    /** Additive methods for use by Evolve. */
    template <typename Base>
    struct EvolveMethods : public Base
    {
      // Note: We don't declare any of the other features until the
      //    derived modules, such as the permittivity/permeability modules.

      /** Terminal case for advancing the H-field. */
      void advanceH(typename Base::TimingIterator it)
      {
        // Nothing to do
      }

      /** Terminal case for writing post H-field advancement. */
      void inspectH(typename Base::TimingIterator it)
      {
        // Nothing to do
      }

      /** Terminal case for advancing the E-field. */
      void advanceE(typename Base::TimingIterator it)
      {
        // Nothing to do
      }

      /** Terminal case for writing post E-field advancement. */
      void inspectE(typename Base::TimingIterator it)
      {
        // Nothing to do
      }
    };
  };

  class FieldBaseBase {};

  //
  // Assemble the field data type
  //

  /** Assembled Field data type. */
  using FieldType = typename Utilities::DecoratedLayers<
      FallbackSelector,
      TimingSelector,
      IndicesSelector,
      FieldDataSelector,
      CoordinateMethodsSelector,
      CoordinatesSelector,
      AdvanceMethodsSelector,
      EvolveMethodsSelector,
      EvolveSelector,
      EvaluateMethodsSelector,
      EvaluateSelector
      >::template For<
          FieldBase, Timing, Indices,
          FieldData, Decorators...>::template Type<FieldBaseBase>;

  /** Instantiation of the internal simulation data. */
  FieldType m_Field;

public:

  /** Method used to run the current instance. */
  void evaluate(void)
  {
    m_Field.evaluate();
  }

  /** Get a handle to the internal simulation data. */
  FieldType& data(void)
  {
    return m_Field;
  }

};

} // namespace Fdtd

#endif // #ifndef FDTD_SIMULATION_HPP

