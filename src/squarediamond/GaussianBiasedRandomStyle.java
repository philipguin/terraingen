package squarediamond;
import java.util.List;
import java.util.ListIterator;
import java.util.Random;

import squarediamond.SquareDiamondArray2DPopulator.ISquareDiamondStyle;
import src.Array2D;


public final class GaussianBiasedRandomStyle extends SimpleCompositeSquareDiamondStyle
{
	private final Random random;
	private final List<Float> randomnesses;
	
	private ListIterator<Float> it;
	private float magnitude;
	
	public GaussianBiasedRandomStyle(Random random, ISquareDiamondStyle base, List<Float> randomnesses)
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
			magnitude = it.next();
	}
	
	@Override
	public final float bias(float value)
	{
		return Math.max(0f, Math.min(1f, value + magnitude * (float)random.nextGaussian()));
	}
}
