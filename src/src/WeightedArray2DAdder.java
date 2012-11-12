package src;

import java.util.List;
import java.util.ListIterator;

public class WeightedArray2DAdder implements IPopulator<Array2D>
{
	private final List<IPopulator<Array2D>> populators;
	private final List<Float> weights;

	public WeightedArray2DAdder(List<IPopulator<Array2D>> populators, List<Float> weights)
	{
		this(populators, weights, -1f);
	}
	
	public WeightedArray2DAdder(List<IPopulator<Array2D>> populators, List<Float> weights, float newTotalWeight)
	{
		this.populators = populators;
		this.weights = weights;
		
		if (weights.size() != populators.size())
			throw new Error("weights and populators must be of same length!");
		
		if (newTotalWeight > 0f)
		{
			float totalWeight = 0f;
			
			for (float weight : weights)
				totalWeight += weight;
			
			for (int i = 0; i < weights.size(); ++i)
				weights.set(i, weights.get(i) / totalWeight * newTotalWeight);
		}
	}
	
	public void populate(Array2D input)
	{
		Array2D buffer = new Array2D(input.width, input.height);
		
		int populatorIndex = 0;
		
		ListIterator<Float> weightIt = weights.listIterator();
		
		for (IPopulator<Array2D> populator : populators)
		{
			float weight = weightIt.next();
			populator.populate(buffer);

			int i, j;
			
			if (populatorIndex == 0)
			{
				for (j = 0; j < input.height; ++j)
				for (i = 0; i < input.width;  ++i)
				{
					input.set(i, j, buffer.get(i, j) * weight);
				}
			}
			else
			{
				for (j = 0; j < input.height; ++j)
				for (i = 0; i < input.width;  ++i)
				{
					input.set(i, j, input.get(i, j) + buffer.get(i, j) * weight);
				}
			}
			
			++populatorIndex;
		}
	}
}
