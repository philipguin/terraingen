package squarediamond;
import java.util.List;
import java.util.ListIterator;
import java.util.Random;

import src.Array2D;


public final class GaussianRandomlyInterpolatedSquareDiamondStyle extends InterpolatedSquareDiamondStyle
{
	private final List<Float> magnitudes;
	
	private ListIterator<Float> it;
	private float magnitude;
	
	public GaussianRandomlyInterpolatedSquareDiamondStyle(Random random, List<Float> magnitudes)
	{
		super(random);
		
		this.magnitudes = magnitudes;
	}
	
	@Override
	public void setup(Array2D values)
	{
		super.setup(values);
		
		it = magnitudes.listIterator();
	}
	
	@Override
	public void onNextIteration()
	{
		super.onNextIteration();
		
		if (it.hasNext())
			magnitude = it.next() / 2f;
	}

	@Override
	protected final float nextWeight()
	{
		return .5f + magnitude * (float)random.nextGaussian();
	}
}
