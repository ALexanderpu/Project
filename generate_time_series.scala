import scala.io._
import util.Random._
import scala.math._


def outputToCSVFile(strFilename:String, arr:Array[Double]) { scala.tools.nsc.io.File(strFilename).writeAll(arr.map(_.toString()).reduce(_ + "\n" + _) + "\n") }

// ok, now a predator-prey variant where there is an additional stock fed by the main state variables

// Lotka–Volterra equations

def dualPredpreyVariantWithNoise(state:Vector[Double])(dt:Double)(hare1NatalityRate:Double, 
                                                                  perLynx1AnnualHare1ProbOfGettingCaught:Double, 
                                                                  lynx1MortalityRate:Double,
                                                                  perHare1NetLynx1BirthRateIncrement:Double,
                                                                  hare2NatalityRate:Double, 
                                                                  perLynx2AnnualHare2ProbOfGettingCaught:Double, 
                                                                  lynx2MortalityRate:Double,
                                                                  perHare2NetLynx2BirthRateIncrement:Double,
                                                                  perLynx2AnnualHare1ProbOfGettingCaught:Double)
                                                                (noiseCoeffHare1:Double, noiseCoeffLynx1:Double, noiseCoeffHare2:Double, noiseCoeffLynx2:Double) : Stream[Vector[Double]]  = 
    {
     // in short, we have
    //  w <==> x
    //  y <==> z
    //  z ==> w  This is the link between the two subsystems.
    def dualPredpreyVariantWithNoiseInternal(state:Vector[Double]) : Stream[Vector[Double]] = state #:: 
        {
        val Vector(hares1,lynx1,hares2,lynx2) = state
        dualPredpreyVariantWithNoiseInternal(Vector(max(0, hares1 + dt * (hare1NatalityRate*hares1 - (perLynx1AnnualHare1ProbOfGettingCaught*lynx1 + perLynx2AnnualHare1ProbOfGettingCaught*lynx2)*hares1 + noiseCoeffHare1 * nextGaussian)),            // x depends on y (and itself)
                                                    max(0, lynx1 + dt * (-lynx1MortalityRate*lynx1 + (perHare1NetLynx1BirthRateIncrement*hares1)*lynx1 + noiseCoeffLynx1 * nextGaussian)),          // y depends on x (and itself)
                                                    max(0, hares2 + dt * (hare2NatalityRate*hares2 - (perLynx2AnnualHare2ProbOfGettingCaught*lynx2)*hares2 + noiseCoeffHare2 * nextGaussian)),            // x depends on y (and itself)
                                                    max(0, lynx2 + dt * (-lynx2MortalityRate*lynx2 + (perHare2NetLynx2BirthRateIncrement*hares2)*lynx2   + noiseCoeffLynx2 * nextGaussian))))
        }
    dualPredpreyVariantWithNoiseInternal(state)
    }


// Now a VERY NOISY  (5x as much as variant 3) case WITH coupling between the two systems  (equal to variant 3) 
// These parameters are established the same as in variant 3, namely to put the natural cycling on different frequencies, with the 
// second system (hare2, lynx2) cycling 1.5-2 times as fast as the first system
// However, there is high noise

// Now a variant with stronger coupling, and HIGH noise
// These parameters are established to put the natural cycling on different frequencies, with the 
// second system (hare2, lynx2) cycling 1.5-2 times as fast as the first system

// build multiple model at least 1000 models
def create_models(num:Int){

  // change initial stocks
  val meanState = Vector(2.0, 4.0, 4.0, 3.0)

  // change parameters: birth rate for prey; death rate for predator
  val meanHare1NatalityRate = 0.02
  val meanHare2NatalityRate = 0.025

  val meanLynx1MortalityRate = 0.03
  val meanLynx2MortalityRate = 0.03

  val meanPerLynx1AnnualHare1ProbOfGettingCaught = 0.005
  val meanPerLynx2AnnualHare2ProbOfGettingCaught = 0.010
  val meanPerLynx2AnnualHare1ProbOfGettingCaught = 0.0025

  val meanPerHare1NetLynx1BirthRateIncrement = 0.0005
  val meanPerHare2NetLynx2BirthRateIncrement = 0.001

  val arraySize = 2000

    //  w <==> x
    //  y <==> z
    //  z  ==> w  
  var item = 0
  for(item <- 1 to num){
    println("This is model " + item)
    val startTime = System.nanoTime
    val state = Vector(max(1, meanState(0) + nextGaussian), max(3, meanState(1) + nextGaussian), max(3, meanState(2) + nextGaussian), max(2, meanState(3) + nextGaussian))
    
    val hare1NatalityRate = meanHare1NatalityRate + 0.01*nextGaussian
    val hare2NatalityRate = meanHare2NatalityRate + 0.01*nextGaussian

    val lynx1MortalityRate = meanLynx1MortalityRate + 0.01*nextGaussian
    val lynx2MortalityRate = meanLynx2MortalityRate + 0.01*nextGaussian

    val perLynx1AnnualHare1ProbOfGettingCaught = meanPerLynx1AnnualHare1ProbOfGettingCaught + 0.002 * nextGaussian
    val perLynx2AnnualHare2ProbOfGettingCaught = meanPerLynx2AnnualHare2ProbOfGettingCaught + 0.003 * nextGaussian
    val perLynx2AnnualHare1ProbOfGettingCaught = meanPerLynx2AnnualHare2ProbOfGettingCaught + 0.001 * nextGaussian

    val perHare1NetLynx1BirthRateIncrement = meanPerHare1NetLynx1BirthRateIncrement + 0.0002 * nextGaussian
    val perHare2NetLynx2BirthRateIncrement = meanPerHare2NetLynx2BirthRateIncrement + 0.0002 * nextGaussian

    val  dualPredpreyVariantStream = dualPredpreyVariantWithNoise(state)(1E-2)(hare1NatalityRate, perLynx1AnnualHare1ProbOfGettingCaught, lynx1MortalityRate, perHare1NetLynx1BirthRateIncrement, hare2NatalityRate, perLynx2AnnualHare2ProbOfGettingCaught, lynx2MortalityRate, perHare2NetLynx2BirthRateIncrement,perLynx2AnnualHare1ProbOfGettingCaught
                                                                                            )(5.0, 5.0, 5.0, 5.0)

    val  dualPredpreyVariant = dualPredpreyVariantStream.zipWithIndex.withFilter(p => (p._2 % 1000) == 0).map(_._1).take(arraySize).toArray
    
    val preStr = "./data/dualPredpreyVariant" + item
    outputToCSVFile(preStr + "_prey1.csv", dualPredpreyVariant.map(_(0)).toArray)
    outputToCSVFile(preStr + "_predator1.csv", dualPredpreyVariant.map(_(1)).toArray)
    outputToCSVFile(preStr + "_prey2.csv", dualPredpreyVariant.map(_(2)).toArray)
    outputToCSVFile(preStr + "_predator2.csv", dualPredpreyVariant.map(_(3)).toArray)  
    val endTime = System.nanoTime
    println("duration: " + (endTime - startTime) / 1e9d)
  }
}

// time use: for each model it will cost almost 6 seconds  1000 models ->almost 2 hours
// memory use: 4 files for a model: 41 * 4 kb   for 1000 models it would be 0.16GB

// start with scala -J-Xmx6g generate_time_series.scala
create_models(10)

