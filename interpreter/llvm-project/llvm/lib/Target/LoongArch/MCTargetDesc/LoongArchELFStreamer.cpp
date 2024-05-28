//===-- LoongArchELFStreamer.cpp - LoongArch ELF Target Streamer Methods --===//
//
// Part of the LLVM Project, under the Apache License v2.0 with LLVM Exceptions.
// See https://llvm.org/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//
//
// This file provides LoongArch specific target streamer methods.
//
//===----------------------------------------------------------------------===//

#include "LoongArchELFStreamer.h"
#include "LoongArchAsmBackend.h"
#include "llvm/BinaryFormat/ELF.h"
#include "llvm/MC/MCAssembler.h"
#include "llvm/MC/MCCodeEmitter.h"
#include "llvm/MC/MCObjectWriter.h"

using namespace llvm;

// This part is for ELF object output.
LoongArchTargetELFStreamer::LoongArchTargetELFStreamer(
    MCStreamer &S, const MCSubtargetInfo &STI)
    : LoongArchTargetStreamer(S) {
  // FIXME: select appropriate ABI.
  setTargetABI(STI.getTargetTriple().isArch64Bit() ? LoongArchABI::ABI_LP64D
                                                   : LoongArchABI::ABI_ILP32D);
}

MCELFStreamer &LoongArchTargetELFStreamer::getStreamer() {
  return static_cast<MCELFStreamer &>(Streamer);
}

void LoongArchTargetELFStreamer::finish() {
  LoongArchTargetStreamer::finish();
  MCAssembler &MCA = getStreamer().getAssembler();
  LoongArchABI::ABI ABI = getTargetABI();

  // Figure out the e_flags.
  //
  // Bitness is already represented with the EI_CLASS byte in the current spec,
  // so here we only record the base ABI modifier. Also set the object file ABI
  // version to v1, as upstream LLVM cannot handle the previous stack-machine-
  // based relocs from day one.
  //
  // Refer to LoongArch ELF psABI v2.01 for details.
  unsigned EFlags = MCA.getELFHeaderEFlags();
  EFlags |= ELF::EF_LOONGARCH_OBJABI_V1;
  switch (ABI) {
  case LoongArchABI::ABI_ILP32S:
  case LoongArchABI::ABI_LP64S:
    EFlags |= ELF::EF_LOONGARCH_ABI_SOFT_FLOAT;
    break;
  case LoongArchABI::ABI_ILP32F:
  case LoongArchABI::ABI_LP64F:
    EFlags |= ELF::EF_LOONGARCH_ABI_SINGLE_FLOAT;
    break;
  case LoongArchABI::ABI_ILP32D:
  case LoongArchABI::ABI_LP64D:
    EFlags |= ELF::EF_LOONGARCH_ABI_DOUBLE_FLOAT;
    break;
  case LoongArchABI::ABI_Unknown:
    llvm_unreachable("Improperly initialized target ABI");
  }
  MCA.setELFHeaderEFlags(EFlags);
}

namespace {
class LoongArchELFStreamer : public MCELFStreamer {
public:
  LoongArchELFStreamer(MCContext &C, std::unique_ptr<MCAsmBackend> MAB,
                       std::unique_ptr<MCObjectWriter> MOW,
                       std::unique_ptr<MCCodeEmitter> MCE)
      : MCELFStreamer(C, std::move(MAB), std::move(MOW), std::move(MCE)) {}
};
} // end namespace

namespace llvm {
MCELFStreamer *createLoongArchELFStreamer(MCContext &C,
                                          std::unique_ptr<MCAsmBackend> MAB,
                                          std::unique_ptr<MCObjectWriter> MOW,
                                          std::unique_ptr<MCCodeEmitter> MCE,
                                          bool RelaxAll) {
  LoongArchELFStreamer *S = new LoongArchELFStreamer(
      C, std::move(MAB), std::move(MOW), std::move(MCE));
  S->getAssembler().setRelaxAll(RelaxAll);
  return S;
}
} // end namespace llvm
