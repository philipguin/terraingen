package squarediamond;
import java.util.List;
import java.util.ListIterator;
import java.util.Random;

import squarediamond.SquareDiamondArray2DPopulator.ISquareDiamondStyle;
import src.Array2D;


public final class UniformWeightedRandomStyle extends SimpleCompositeSquareDiamondStyle
{
	private final Random random;
	private final List<Float> randomnesses;
	
	private ListIterator<Float> it;
	private float randomness, inverseRandomness;
	
	public UniformWeightedRandomStyle(Random random, ISquareDiamondStyle base, List<Float> randomnesses)
	{
		super(base);
		
		this.random = random;
		this.randomnesses = randomnesses;
	}
	
	@Override
	public void setup(Array2D values)
	{
		super.setup(values);
		
		it = randomnesses.listIterator();
	}
	
	@Override
	public void onNextIteration()
	{
		super.onNextIteration();
		
		if (it.hasNext())
		{
			randomness = it.next();
			inverseRandomness = 1f - randomness;
		}
	}
	
	@Override
	public final float bias(float value)
	{
		return random.nextFloat() * randomness + inverseRandomness * value;
	}
}
