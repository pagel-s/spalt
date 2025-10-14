#include <gtest/gtest.h>
#include <filesystem>
#include <fstream>
#include <string>

#include "molecule.h"
#include "surface.h"

TEST(Surface, GenerateSurface_OpensFile) {
    const auto tmpPath = std::filesystem::temp_directory_path() / "spalt_test.xyzr";

    {
        std::ofstream out(tmpPath);
        ASSERT_TRUE(out.good()) << "Failed to create temp xyzr file: " << tmpPath;
        out << R"(  27.331   -1.048   34.864 2.0000 1 A_1_x_PRO_N_x
  27.842    0.076   34.073 2.0000 1 A_1_x_PRO_CA_x
  26.776    1.124   33.768 1.7000 1 A_1_x_PRO_C_x
  25.604    0.783   33.611 1.4000 1 A_1_x_PRO_O_x
  28.328   -0.592   32.775 2.0000 1 A_1_x_PRO_CB_x
  27.950   -2.062   32.878 2.0000 1 A_1_x_PRO_CG_x
  26.992   -2.198   34.013 2.0000 1 A_1_x_PRO_CD_x
  28.639    0.503   34.536 1.0000 1 A_1_x_PRO_HA_x
  27.886   -0.175   31.972 1.0000 1 A_1_x_PRO_HB2_x
  29.325   -0.495   32.675 1.0000 1 A_1_x_PRO_HB3_x
  27.536   -2.366   32.006 1.0000 1 A_1_x_PRO_HG2_x
  28.784   -2.614   33.029 1.0000 1 A_1_x_PRO_HG3_x
  26.042   -2.138   33.707 1.0000 1 A_1_x_PRO_HD2_x
  27.127   -3.055   34.513 1.0000 1 A_1_x_PRO_HD3_x
  26.399   -0.853   35.164 1.0000 1 A_1_x_PRO_H_x
  27.918   -1.196   35.656 1.0000 1 A_1_x_PRO_H2_x
  27.190    2.390   33.691 1.5000 1 A_2_x_THR_N_x
  26.297    3.471   33.296 2.0000 1 A_2_x_THR_CA_x
  26.498    3.904   31.852 1.7000 1 A_2_x_THR_C_x
  25.605    4.541   31.282 1.4000 1 A_2_x_THR_O_x
  26.486    4.683   34.217 2.0000 1 A_2_x_THR_CB_x
  27.848    5.124   34.164 1.4000 1 A_2_x_THR_OG1_x
  26.127    4.321   35.651 2.0000 1 A_2_x_THR_CG2_x
  28.191    2.532   33.935 1.0000 1 A_2_x_THR_H_x
  25.340    3.149   33.401 1.0000 1 A_2_x_THR_HA_x
  25.893    5.446   33.887 1.0000 1 A_2_x_THR_HB_x
  28.350    4.515   33.522 1.0000 1 A_2_x_THR_HG1_x
  26.403    5.062   36.268 1.0000 1 A_2_x_THR_HG21_x
  25.137    4.180   35.729 1.0000 1 A_2_x_THR_HG22_x
  26.600    3.478   35.917 1.0000 1 A_2_x_THR_HG23_x
)";
    }

    std::cout << "tmpPath: " << tmpPath << std::endl;

    Molecule molecule("C1CCCCC1");
    // Generate conformers first
    molecule.calculate_conformers(1);
    // Test FPS subsampling by creating surface with fps sampling method
    int conformer_id = molecule.createSurface(0, 100, 1.2, 3.0, 3.0, "tses", "msms", "fps");
    EXPECT_EQ(molecule.getSurfaceVertexCount(conformer_id), 100);

    // Test random subsampling - create a new molecule to avoid caching
    Molecule molecule2("C1CCCCC1OC");
    molecule2.calculate_conformers(1);
    int conformer_id2 = molecule2.createSurface(0, 50, 1.2, 3.0, 3.0, "tses", "msms", "random");
    EXPECT_EQ(molecule2.getSurfaceVertexCount(conformer_id2), 50);

    // Test full (no subsampling) - create a new molecule to avoid caching
    Molecule molecule3("C1CCCCC1");
    molecule3.calculate_conformers(1);
    int conformer_id3 = molecule3.createSurface(0, 1000, 1.2, 3.0, 3.0, "tses", "msms", "full");
    EXPECT_GT(molecule3.getSurfaceVertexCount(conformer_id3), 100);

    std::error_code ec;
    std::filesystem::remove(tmpPath, ec);
    std::filesystem::remove(tmpPath.string() + ".msms", ec);
}
